
#include <stdio.h>
#include <stdlib.h>
#include "bpoly.h"

/* Monomials are stored as a sorted double-linked list. */
typedef struct mon_s mon_t;

struct mon_s{
    bvar_t    x;
    mon_t *prev;
    mon_t *next;
};

struct bpoly_s{
    mon_t *mon;
};

bpoly_t *bpoly_new_zero()
{
    bpoly_t *bpoly = malloc(sizeof(bpoly_t));
    bpoly->mon = NULL;
    
    return bpoly;
}

bpoly_t *bpoly_new_one()
{    
    return bpoly_new((bvar_t[]){0}, 1);
}

bpoly_t *bpoly_new(bvar_t *mon, int num)
{
    bpoly_t *bpoly = malloc(sizeof(bpoly_t));
    bpoly->mon = NULL;
    
    int i;
    for (i = 0; i < num; i++, mon++)
    {
        bpoly_add_mon(bpoly, *mon);
    }
    
    return bpoly;
}

bpoly_t *bpoly_new_random(int n, int t)
{
    bpoly_t *bpoly = malloc(sizeof(bpoly_t));
    bpoly->mon = NULL;
    
    int i;
    for (i = 0; i < t; i++)
    {
        bpoly_add_mon(bpoly, rand() & (((bvar_t)1 << n) - 1));
    }
    
    return bpoly;    
}

void bpoly_free(bpoly_t *bpoly)
{
    mon_t *mon = bpoly->mon;
    
    while (mon)
    {
        mon_t *prev = mon;
        mon = mon->next;
        free(prev);
    }
    
    free(bpoly);
}

bpoly_t *bpoly_copy(bpoly_t *bpoly)
{
    bpoly_t *copy = malloc(sizeof(bpoly_t));
    
    mon_t *mon = bpoly->mon;
    
    if (mon == NULL)
    {
        copy->mon = NULL;
    }
    else
    {
        mon_t *prev = NULL;
        
        while (mon)
        {
            mon_t *mon_copy = malloc(sizeof(mon_t));
            mon_copy->x = mon->x;
            mon_copy->prev = prev;
            
            if (prev)
            {
                prev->next = mon_copy;
            }
            else
            {
                copy->mon = mon_copy;
            }
            
            mon_copy->next = NULL;
            
            prev = mon_copy;
            mon = mon->next;
        }
    }
    
    return copy;
}

void bpoly_print(bpoly_t *bpoly)
{
    mon_t *mon = bpoly->mon;
    
    bool plus = 0;

    if (mon == NULL)
    {
        printf("0");
    }
    else
    {
        while (mon)
        {
            bvar_t x = mon->x;
            
            if (plus)
            {
                printf("+ ");
            }
            
            if (x == 0)
            {
                printf("1 ");
            }
            else
            {
                int i;
                for (i = 0; i < BVAR_MAX; i++)
                {
                    if ((x >> i) & 1)
                    {
                        printf("x_%d ", i + 1);
                    }
                }
            }
            
            plus = 1;
            mon = mon->next;
        }
    }
    
    printf("\n");
}

int bpoly_degree(bpoly_t *bpoly)
{
    if (bpoly->mon == NULL)
    {
        return -1;
    }
    
    mon_t *mon = bpoly->mon;
    int degree = 0;
    
    while (mon)
    {
        int d = __builtin_popcountl(mon->x);
        
        if (d > degree)
        {
            degree = d;
        }
        
        mon = mon->next;
    }
    
    return degree;
    
}

bvar_t bpoly_vars(bpoly_t *bpoly)
{
    bvar_t x = 0;
    mon_t *mon = bpoly->mon;

    while (mon)
    {
        x |= mon->x;
        mon = mon->next;
    }

    return x;
}

bool bpoly_eval(bpoly_t *bpoly, bvar_t x)
{
    bool v = 0;
    mon_t *mon = bpoly->mon;
    
    while (mon)
    {
        v ^= ((mon->x & x) == mon->x);
        mon = mon->next;
    }
    
    return v;
}

void bpoly_add_mon(bpoly_t *bpoly, bvar_t x)
{
    if (bpoly->mon == NULL)
    {
        mon_t *mon = malloc(sizeof(mon_t));
        mon->x = x;
        mon->prev = NULL;
        mon->next = NULL;
        bpoly->mon = mon;
        return;
    }
    
    mon_t *list = bpoly->mon;
    
    while (list->x < x && list->next != NULL)
    {
        list = list->next;
    }
    
    if (list->x == x)
    {
        /* Characteristic 2. */
        if (bpoly->mon == list)
        {
            bpoly->mon = list->next;
        }
        
        if (list->prev)
        {
            list->prev->next= list->next;
        }
        
        if (list->next)
        {
            list->next->prev = list->prev;
        }
        
        free(list);
    }
    else if (list->x > x)
    {
        mon_t *mon = malloc(sizeof(mon_t));
        mon->x = x;
        mon->prev = list->prev;
        mon->next = list;
        if (list->prev)
        {
            list->prev->next = mon;
        }
        list->prev = mon;
        
        if (bpoly->mon == list)
        {
            bpoly->mon = mon;
        }
    }
    else
    {
        /* End of the list. */
        mon_t *mon = malloc(sizeof(mon_t));
        mon->x = x;
        mon->prev = list;
        mon->next = NULL;
        list->next = mon;
    }

}

void bpoly_mul_mon(bpoly_t *bpoly, bvar_t x)
{
    bpoly_t *temp = bpoly_new(NULL, 0);
    
    mon_t *mon = bpoly->mon;
    
    while (mon)
    {
        mon_t *prev = mon;
        bpoly_add_mon(temp, mon->x | x);
        free(prev);
        mon = mon->next;
    }
    
    bpoly->mon = temp->mon;
    free(temp);
}

void bpoly_add(bpoly_t *bpoly1, bpoly_t *bpoly2)
{
    mon_t *mon = bpoly2->mon;
    
    while (mon)
    {
        bpoly_add_mon(bpoly1, mon->x);
        mon = mon->next;
    }
}

void bpoly_mul(bpoly_t *bpoly1, bpoly_t *bpoly2)
{
    bpoly_t *temp = bpoly_new_zero();
    
    mon_t *mon1 = bpoly1->mon;
    
    /* Place product of 'bpoly1' and 'bpoly2' in 'temp' */
    while (mon1)
    {
        mon_t *mon2 = bpoly2->mon;
        
        while (mon2)
        {
            bpoly_add_mon(temp, mon1->x | mon2->x);
            mon2 = mon2->next;
        }
        
        mon1 = mon1->next;
    }

    /* Destroy monomials of 'bpoly1' */
    mon1 = bpoly1->mon;
    while (mon1)
    {
        mon_t *prev = mon1;
        free(prev);
        mon1 = mon1->next;
    }
    
    /* Replace 'bpoly1' with 'temp' */
    bpoly1->mon = temp->mon;
    free(temp);
}

void bpoly_assign(bpoly_t *bpoly, bvar_t x, bvar_t a)
{
    bpoly_t *temp = bpoly_new(NULL, 0);
    
    mon_t *mon = bpoly->mon;
    
    while (mon)
    {    
        bvar_t y = mon->x & x;
        if ((a & y) == y){
            bpoly_add_mon(temp, mon->x & ~x);
        }
        
        mon = mon->next;
    }
    
    bpoly->mon = temp->mon;
    free(temp);
}

bfunc_t *bpoly_monomials(bpoly_t *bpoly, int n)
{
    bfunc_t *bfunc = bfunc_new(n);
    
    mon_t *mon = bpoly->mon;
    
    while (mon)
    {
        bfunc_set(bfunc, mon->x, 1);
        mon = mon->next;
    }

    return bfunc;
}

bfunc_t *bpoly_truth_table(bpoly_t *bpoly, int n)
{
    bfunc_t *bfunc = bpoly_monomials(bpoly, n);
    
    bfunc_zeta_transform(bfunc);
    
    return bfunc;
}
