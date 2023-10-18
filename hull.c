#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define BUFFSIZE 100

struct vec {
    double x;
    double y;
};

struct vecset {
    struct vec *data;
    size_t size;
    size_t capacity;
};

const int A = 2;

// Functions on vectors
/**
 * Allocates memory and creates a vector
 * @param x The x coordinate
 * @param y The y coordinate
 * @return A pointer to the newly created vector
 */
struct vec *vec_create(double x, double y) {
    struct vec *new = malloc(sizeof(struct vec));
    new->x = x;
    new->y = y;
    return new;
}

/**
 * Returns the dot of two vectors
 * @param v1 The first vector
 * @param v2 The second vector
 * @return The dot of the two vectors
 */
double dot(struct vec *v1, struct vec *v2) {
    return v1->x * v2->x + v1->y * v2->y;
}

/**
 * Returns the cross product between three vectors
 * @param v1 The first vector
 * @param v2 The second vector
 * @param v3 The third vector
 * @return The cross product of the vectors
 */
double cross(struct vec *v1, struct vec *v2, struct vec *v3) {
    return (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
}

/**
 * Checks if three vectors are in a left turn (cross product strictly negative)
 * @param v1 The first vector
 * @param v2 The second vector
 * @param v3 The third vector
 * @return true if the vectors are in a left turn, false otherwise
 */
bool is_left_turn(struct vec *v1, struct vec *v2, struct vec *v3) {
    return cross(v1, v2, v3) < 0;
}

// Functions on vector sets
/**
 * Creates a vecset and allocated its default capacity
 * @param self The vecset to create
 */
void vecset_create(struct vecset *self) {
    const int defaultCapacity = 10;
    self->capacity = defaultCapacity;
    self->size = 0;
    self->data = calloc(defaultCapacity, sizeof(struct vec));
}

/**
 * Frees the memory of a vecset
 * @param self The vecset to destroy
 */
void vecset_destroy(struct vecset *self) {
    free(self->data);
}

/**
 * Allocated more memory and increases the capacity of a vecset
 * @param self The vecset to grow
 */
void vecset_grow(struct vecset *self) {
    self->capacity *= A;
    struct vec *new = calloc(self->capacity, sizeof(struct vec));
    memcpy(new, self->data, self->size * sizeof(struct vec));
    free(self->data);
    self->data = new;
}

/**
 * Adds a vector to a vecset
 * @param self The vecset to add the vector in
 * @param p The vector to add to the vecset
 */
void vecset_add(struct vecset *self, struct vec *p) {
    if (self->size == self->capacity) vecset_grow(self);
    self->data[self->size] = *p;
    self->size++;
}

typedef int (*comp_func_t) (const struct vec *p1, const struct vec *p2, const void *ctx);

/**
 * Returns the maximum vector according to a comparison function
 * @param self The vecset to search in
 * @param func The comparison function
 * @param ctx The context for the comparison function
 * @return A pointer to the maximum vector
 */
const struct vec *vecset_max(const struct vecset *self, comp_func_t func, const void *ctx) {
    assert(self->size > 0);
    struct vec *max = &self->data[0];
    for (size_t i = 0; i < self->size; i++) {
        if (func(max, &self->data[i], ctx) > 0) max = &self->data[i];
    }
    return max;
}

/**
 * Returns the minimum vector according to a comparison function
 * @param self The vecset to search in
 * @param func The comparison function
 * @param ctx The context for the comparison function
 * @return A pointer to the minimum vector
 */
const struct vec *vecset_min(const struct vecset *self, comp_func_t func, const void *ctx) {
    assert(self->size > 0);
    struct vec *min = &self->data[0];
    for (size_t i = 0; i < self->size; i++) {
        if (func(min, &self->data[i], ctx) < 0) min = &self->data[i];
    }
    return min;
}

/**
 * Merges two sorted arrays of vectors in one sorted array.
 * The two original arrays are freed.
 * The result array must be already allocated, and is assumed to have enough memory to store both functions (size equal
 * to size1 + size2)
 * @param a1 The first sorted array
 * @param size1 The size of the first array
 * @param a2 The second sorted array
 * @param size2 The size of the second array
 * @param func A comparison function used to compare vectors
 * @param ctx The context for the comparison function
 */
void vec_array_merge(struct vec *a1, size_t size1, struct vec *a2, size_t size2, struct vec *result,
        comp_func_t func, const void *ctx) {
    size_t index = 0;
    while (size1 > 0 || size2 > 0) {
        if (size2 == 0 || func(&a1[size1 - 1], &a2[size2 - 1], ctx) < 0) {
            result[index] = a1[size1 - 1];
            size1--;
        }
        else {
            result[index] = a2[size2 - 1];
            size2--;
        }
        index++;
    }
    free(a1);
    free(a2);
}

/**
 * Splits an array of vectors in two smaller arrays.
 * The smaller arrays must be already allocated, and are assumed to have enough memory to store the result.
 * (size1 = size / 2 and size2 = size - size / 2)
 * @param a1 The array storing the first half of the original array
 * @param a2 The array storing the second half of the original array
 * @param split The array to split
 * @param size The size of the array
 */
void vec_array_split(struct vec *a1, struct vec *a2, struct vec *split, size_t size) {
    size_t index1 = 0, index2 = 0;
    size_t i = 0;
    for (; i < size / 2; i++) {
        a1[index1] = split[i];
        index1++;
    }
    for (; i < size; i++) {
        a2[index2] = split[i];
        index2++;
    }
}

/**
 * Sorts an array of vectors using a merge sort, according to a comparison function.
 * @param array The array to sort
 * @param size The size of the array
 * @param func A comparison function
 * @param ctx The context for the comparison function
 */
void vec_array_merge_sort(struct vec *array, size_t size, comp_func_t func, const void *ctx) {
    if (size == 1) return;
    if (size == 2) {
        if (func(&array[0], &array[1], ctx) > 0) {
            struct vec temp = array[0];
            array[0] = array[1];
            array[1] = temp;
        }
        return;
    }
    size_t size1 = size - size / 2;
    size_t size2 = size / 2;
    struct vec *a1 = calloc(size1, sizeof(struct vec));
    struct vec *a2 = calloc(size2, sizeof(struct vec));
    vec_array_split(a1, a2, array, size);
    vec_array_merge(a1, size1, a2, size2, array, func, ctx);
}

/**
 * Sorts a vecset using a merge sort, according to a comparison function
 * @param self The vecset to sort
 * @param func The comparison function
 * @param ctx The context for the comparison function
 */
void vecset_sort(struct vecset *self, comp_func_t func, const void *ctx) {
    vec_array_merge_sort(self->data, self->size, func, ctx);
}

/**
 * Adds a vector at the end of the list in a vecset
 * @param self The vecset to add the vector in
 * @param p The vector to add in the vecset
 */
void vecset_push(struct vecset *self, struct vec p) {
    vecset_add(self, &p);
}

void vecset_pop(struct vecset *self) {
    assert(self->size > 0);
    self->size--;
}

const struct vec *vecset_top(const struct vecset *self) {
    assert(self->size > 0);
    return &self->data[self->size - 1];
}

const struct vec *vecset_second(const struct vecset *self) {
    assert(self->size > 1);
    return &self->data[self->size - 2];
}

int main() {
    setbuf(stdout, NULL);

    char buffer[BUFFSIZE];
    fgets(buffer, BUFFSIZE, stdin);

    size_t count = strtol(buffer, NULL, 10);

    for (size_t i = 0; i < count; ++i) {
        struct vec p;

        fgets(buffer, BUFFSIZE, stdin);

        char *endptr = buffer;
        p.x = strtod(endptr, &endptr);
        p.y = strtod(endptr, &endptr);
    }
    return 0;
}
