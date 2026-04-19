/* elements.h — static element symbol table (Z = 1..118).
 *
 * Provides symbol-to-atomic-number and atomic-number-to-symbol lookups.
 * Include this header in exactly one translation unit that needs it.
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <string.h>
#include <ctype.h>

#define NUM_ELEMENTS 119  /* index 0 unused; 1..118 */

static const char *ELEMENT_SYMBOLS[NUM_ELEMENTS] = {
    "",                                                         /*  0 */
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", /*  1-10  */
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", /* 11-20  */
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", /* 21-30  */
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", /* 31-40  */
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", /* 41-50  */
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", /* 51-60  */
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", /* 61-70  */
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", /* 71-80  */
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", /* 81-90  */
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", /* 91-100 */
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", /*101-110 */
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"              /*111-118 */
};

/* Case-insensitive symbol lookup.  Returns atomic number (1..118) or -1.
 * Marked unused because elements.h is included by translation units
 * (e.g. converter.c) that only need ELEMENT_SYMBOLS[], not this helper. */
__attribute__((unused))
static int element_by_symbol(const char *sym)
{
    char buf[4] = {0};
    int  len = 0;
    while (sym[len] && len < 3) {
        buf[len] = (len == 0) ? (char)toupper((unsigned char)sym[len])
                              : (char)tolower((unsigned char)sym[len]);
        len++;
    }
    buf[len] = 0;
    for (int i = 1; i < NUM_ELEMENTS; i++)
        if (strcmp(buf, ELEMENT_SYMBOLS[i]) == 0) return i;
    return -1;
}

#endif /* ELEMENTS_H */
