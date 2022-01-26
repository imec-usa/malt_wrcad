// vi: ts=2 sts=2 sw=2 et tw=100

#include "list.h"
#include <stdlib.h>

list_t EMPTY_LIST = {NULL, 0, 0};

/* Reverses the order of strings in a list. */
void lst_reverse(list_t *lst)
{
  for (size_t l = 0, r = lst->len; l < --r; ++l) {
    char *temp = lst->ptr[l];
    lst->ptr[l] = lst->ptr[r];
    lst->ptr[r] = temp;
  }
}

/* Adds a new string to the end of a list, reallocating if necessary. */
list_t *lst_push(list_t *lst, char *str)
{
  size_t end = lst->len++;
  if (lst->len > lst->cap) {
    // growth sequence: 0, 4, 10, 19, 32, 52, ...
    lst->cap = lst->cap * 3 / 2 + 4;
    lst->ptr = realloc(lst->ptr, lst->cap * sizeof *lst->ptr);
  }
  lst->ptr[end] = str;
  return lst;
}

/* Frees `*lst` and all its contents. */
void lst_drop(list_t *lst)
{
  for (size_t i = 0; i < lst->len; ++i) {
    free(lst->ptr[i]);
  }
  free(lst->ptr);
  *lst = EMPTY_LIST;
}

bool lst_empty(list_t *lst) { return lst->len == 0; }

/* Returns the last element of `*lst`. */
char *lst_last(list_t *lst) { return lst->ptr[lst->len - 1]; }
