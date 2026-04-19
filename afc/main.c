/* main.c — CLI for the Atomic Format Converter (afc).
 *
 * Usage: afc [options] input output
 *
 *   -i FMT   Input format  (xyz, lammps, cel).  Auto-detected if omitted.
 *   -o FMT   Output format (xyz, lammps, cel).  Auto-detected if omitted.
 *   -t MAP   Type map, e.g. "1=Ni,2=Al" or "Ni,Al".
 *   -h       Help.
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#include "converter.h"

#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options] input output\n"
        "\n"
        "Convert between Extended XYZ, LAMMPS data (atomic), and CEL formats.\n"
        "\n"
        "Options:\n"
        "  -i FMT   Input format:  xyz, lammps, cel  (auto-detected by extension)\n"
        "  -o FMT   Output format: xyz, lammps, cel  (auto-detected by extension)\n"
        "  -t MAP   Type map, e.g. \"1=Ni,2=Al\" or \"Ni,Al\" (types assigned in order)\n"
        "  -h       Show this help\n"
        "\n"
        "Supported extensions:\n"
        "  .xyz[.gz]                  Extended XYZ\n"
        "  .lmp[.gz] .data[.gz]       LAMMPS data (atomic style)\n"
        "  .cel[.gz]                  CEL (Dr Probe supercell)\n"
        "\n"
        "The type map (-t) is required when converting to CEL from formats\n"
        "that use numeric type IDs. Example: -t Ni or -t 1=Ni,2=Al\n",
        prog);
}

/* Parse a type map string.  Accepted forms:
 *   "Ni"            →  type 1 = Ni
 *   "Ni,Al"         →  type 1 = Ni, type 2 = Al
 *   "1=Ni,2=Al"     →  explicit mapping
 */
static int parse_typemap(const char *str, TypeMap *tm)
{
    const char *p = str;
    int next_type = 1;

    while (*p) {
        /* skip commas and spaces */
        while (*p == ',' || *p == ' ') p++;
        if (!*p) break;

        int type_id = next_type;

        /* check for explicit "N=" prefix */
        const char *eq = strchr(p, '=');
        const char *comma = strchr(p, ',');
        if (eq && (!comma || eq < comma)) {
            type_id = atoi(p);
            if (type_id < 1 || type_id >= AFC_MAX_TYPES) {
                fprintf(stderr, "afc: invalid type id %d in type map\n", type_id);
                return -1;
            }
            p = eq + 1;
        }

        /* read symbol */
        char sym[AFC_MAX_SYMBOL] = {0};
        int i = 0;
        while (*p && *p != ',' && *p != ' ' && i < AFC_MAX_SYMBOL - 1)
            sym[i++] = *p++;
        sym[i] = 0;

        if (!sym[0]) {
            fprintf(stderr, "afc: empty symbol in type map\n");
            return -1;
        }

        if (type_id >= AFC_MAX_TYPES) {
            fprintf(stderr, "afc: type %d exceeds maximum\n", type_id);
            return -1;
        }
        strcpy(tm->symbols[type_id], sym);
        if (type_id >= tm->ntypes) tm->ntypes = type_id;
        next_type = type_id + 1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    setlocale(LC_NUMERIC, "C");

    const char *in_fmt_str  = NULL;
    const char *out_fmt_str = NULL;
    const char *typemap_str = NULL;
    const char *input_path  = NULL;
    const char *output_path = NULL;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            in_fmt_str = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            out_fmt_str = argv[++i];
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            typemap_str = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "afc: unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        } else if (!input_path) {
            input_path = argv[i];
        } else if (!output_path) {
            output_path = argv[i];
        } else {
            fprintf(stderr, "afc: unexpected argument: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (!input_path || !output_path) {
        fprintf(stderr, "afc: input and output files are required\n\n");
        usage(argv[0]);
        return 1;
    }

    /* Determine formats */
    Format in_fmt = in_fmt_str ? format_from_string(in_fmt_str)
                               : format_from_extension(input_path);
    Format out_fmt = out_fmt_str ? format_from_string(out_fmt_str)
                                 : format_from_extension(output_path);

    if (in_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "afc: cannot determine input format for %s. "
                "Use -i to specify.\n", input_path);
        return 1;
    }
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "afc: cannot determine output format for %s. "
                "Use -o to specify.\n", output_path);
        return 1;
    }

    /* Read */
    Config cfg = config_create();

    fprintf(stderr, "afc: reading %s (%s) ...\n", input_path, format_name(in_fmt));
    if (config_read(&cfg, input_path, in_fmt) != 0) {
        fprintf(stderr, "afc: failed to read %s\n", input_path);
        config_free(&cfg);
        return 1;
    }
    fprintf(stderr, "afc: %d atoms, %d types\n",
            cfg.atoms.count, cfg.typemap.ntypes);

    /* Apply user type map (overrides anything from file) */
    if (typemap_str) {
        TypeMap user_tm;
        memset(&user_tm, 0, sizeof(user_tm));
        if (parse_typemap(typemap_str, &user_tm) != 0) {
            config_free(&cfg);
            return 1;
        }
        /* Merge: user-specified entries override file-derived ones */
        for (int t = 1; t <= user_tm.ntypes; t++) {
            if (user_tm.symbols[t][0])
                strcpy(cfg.typemap.symbols[t], user_tm.symbols[t]);
        }
        if (user_tm.ntypes > cfg.typemap.ntypes)
            cfg.typemap.ntypes = user_tm.ntypes;
    }

    /* Write */
    fprintf(stderr, "afc: writing %s (%s) ...\n", output_path, format_name(out_fmt));
    if (config_write(&cfg, output_path, out_fmt) != 0) {
        fprintf(stderr, "afc: failed to write %s\n", output_path);
        config_free(&cfg);
        return 1;
    }

    fprintf(stderr, "afc: done.\n");
    config_free(&cfg);
    return 0;
}
