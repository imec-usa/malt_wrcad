// vi: ts=3 sts=2 sw=2 et tw=100

#pragma once

#include <stdbool.h>
#include <stddef.h>  // for size_t

struct list {
  char **ptr;
  size_t len;
  size_t cap;
};

typedef struct list list_t;

extern list_t EMPTY_LIST;

void lst_reverse(list_t *lst);
list_t *lst_push(list_t *lst, char *str);
char *lst_pop(list_t *lst);
void lst_drop(list_t *lst);
bool lst_empty(list_t *lst);
char *lst_last(list_t *lst);
