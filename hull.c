#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define BUFFSIZE 100

/**
 * A 2D vector
 */
struct vec {
    double x;
    double y;
};

/**
 * A structure representing a cloud of points as a dynamic list
 */
struct vecset {
    struct vec *data;
    size_t size;
    size_t capacity;
};

const int A = 2;

// Functions on vectors
/**
 * Creates a vector
 * @param x The x coordinate
 * @param y The y coordinate
 * @return The newly created vector
 */
struct vec vec_create(double x, double y) {
    struct vec new;
    new.x = x;
    new.y = y;
    return new;
}

/**
 * Prints a vector in the console
 * @param self The vector to print
 */
void vec_print(struct vec *self) {
    printf("x: %f, y: %f", self->x, self->y);
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
    return cross(v1, v2, v3) > 0;
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

void vecset_print(struct vecset *self) {
    for (size_t i = 0; i < self->size; i++) {
        vec_print(&self->data[i]);
        printf("\n");
    }
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
void vecset_add(struct vecset *self, struct vec p) {
    if (self->size == self->capacity) vecset_grow(self);
    self->data[self->size] = p;
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
    for (size_t i = 1; i < self->size; i++) {
        if (func(max, &self->data[i], ctx) < 0) max = &self->data[i];
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
    for (size_t i = 1; i < self->size; i++) {
        if (func(min, &self->data[i], ctx) > 0) min = &self->data[i];
    }
    return min;
}

/**
 * Merges two sorted arrays of vectors in one sorted array.
 * The two original arrays are freed.
 * The result array must be already allocated, and is assumed to have enough memory to store the result (size equal
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
    size_t index1 = 0;
    size_t index2 = 0;
    while (index1 < size1 || index2 < size2) {
        if (index1 < size1 && (index2 == size2 || func(&a1[index1], &a2[index2], ctx) < 0)) {
            result[index] = a1[index1++];
        }
        else result[index] = a2[index2++];
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
    size_t size1 = size / 2;
    size_t size2 = size / 2 + size % 2;
    struct vec *a1 = calloc(size1, sizeof(struct vec));
    struct vec *a2 = calloc(size2, sizeof(struct vec));
    vec_array_split(a1, a2, array, size);
    vec_array_merge_sort(a1, size1, func, ctx);
    vec_array_merge_sort(a2, size2, func, ctx);
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
    vecset_add(self, p);
}

/**
 * Returns the second element of a vecset from the end
 * @param self The vecset to get the vector from
 * @return The second last vector of the vecset
 */
const struct vec *vecset_second(const struct vecset *self) {
    assert(self->size > 1);
    return &self->data[self->size - 2];
}

/**
 * Removes the last vector of a vecset
 * @param self The vecset to remove a vector from
 */
void vecset_pop(struct vecset *self) {
    assert(self->size > 0);
    self->size--;
}

/**
 * Returns the element at the end of a vecset
 * @param self The vecset to get the vector from
 * @return The last vector of the vecset
 */
const struct vec *vecset_top(const struct vecset *self) {
    assert(self->size > 0);
    return &self->data[self->size - 1];
}

int comp_distance_to_origin(const struct vec *p1, const struct vec *p2, const void *ctx) {
    return (int)(p1->x + p1->y - p2->x - p2->y);
}

void jarvis_march(const struct vecset *in, struct vecset *out) {
}

void tests() {
    // Test vec
    struct vec vec = vec_create(10, 20);
    assert(vec.x == 10);
    assert(vec.y == 20);

    // Test dot
    struct vec vec1 = vec_create(10, 20);
    struct vec vec2 = vec_create(40, 80);
    assert(dot(&vec1, &vec2) == (10 * 40 + 20 * 80));

    // Test cross
    vec1 = vec_create(10, 60);
    vec2 = vec_create(978, 68);
    struct vec vec3 = vec_create(753, 354);
    assert(cross(&vec1, &vec2, &vec3) == (978 - 10) * (354 - 60) - (68 - 60) * (753 - 10));

    // Test left turn
    vec1 = vec_create(0, 0);
    vec2 = vec_create(0, 10);
    vec3 = vec_create(10, 5);
    assert(!is_left_turn(&vec1, &vec2, &vec3));
    vec3 = vec_create(-10, 5);
    assert(is_left_turn(&vec1, &vec2, &vec3));

    // Test vecset
    struct vecset vecset;
    vecset_create(&vecset);
    assert(vecset.capacity == 10);
    assert(vecset.size == 0);
    assert(vecset.data != NULL);

    // Test vecset_grow
    vecset_grow(&vecset);
    assert(vecset.capacity == 20);
    vecset_destroy(&vecset);

    // Test vecset_add
    vecset_create(&vecset);
    struct vec vec4 = vec_create(10, 10);
    vecset_add(&vecset, vec4);
    assert(vecset.data[0].x == 10);
    assert(vecset.data[0].y == 10);
    assert(vecset.size == 1);
    for (int i = 0; i < 2000000; i++) {
        struct vec vec5 = vec_create(i, i);
        vecset_add(&vecset, vec5);
        assert(vecset.size == 2 + i);
    }
    vecset_destroy(&vecset);

    // Test vecset_max
    vec1 = vec_create(0, 0);
    vec2 = vec_create(5, 5);
    vec3 = vec_create(10, 10);
    struct vecset vecset1;
    vecset_create(&vecset1);
    vecset_add(&vecset1, vec1);
    vecset_add(&vecset1, vec2);
    vecset_add(&vecset1, vec3);
    const struct vec *max = vecset_max(&vecset1, &comp_distance_to_origin, NULL);
    assert(max->x == 10 && max->y == 10);

    // Test vecset_min
    const struct vec *min = vecset_min(&vecset1, &comp_distance_to_origin, NULL);
    assert(min->x == 0 && min->y == 0);

    // Test vecset_sort
    // Already sorted
    printf("Already sorted:\n");
    vecset_sort(&vecset1, &comp_distance_to_origin, NULL);
    printf("After sort:\n");
    vecset_print(&vecset1);
    printf("\n");
    assert(vecset1.data[0].x == 0 && vecset1.data[0].y == 0);
    assert(vecset1.data[1].x == 5 && vecset1.data[1].y == 5);
    assert(vecset1.data[2].x == 10 && vecset1.data[2].y == 10);
    vecset_destroy(&vecset1);
    // Inverted
    vecset_create(&vecset1);
    vecset_add(&vecset1, vec3);
    vecset_add(&vecset1, vec2);
    vecset_add(&vecset1, vec1);
    printf("Inverted:\n");
    vecset_sort(&vecset1, &comp_distance_to_origin, NULL);
    printf("After sort:\n");
    vecset_print(&vecset1);
    printf("\n");
    assert(vecset1.data[0].x == 0 && vecset1.data[0].y == 0);
    assert(vecset1.data[1].x == 5 && vecset1.data[1].y == 5);
    assert(vecset1.data[2].x == 10 && vecset1.data[2].y == 10);
    vecset_destroy(&vecset1);
    // Random order
    vecset_create(&vecset1);
    vecset_add(&vecset1, vec3);
    vecset_add(&vecset1, vec1);
    vecset_add(&vecset1, vec2);
    printf("Random order:\n");
    vecset_sort(&vecset1, &comp_distance_to_origin, NULL);
    printf("After sort:\n");
    vecset_print(&vecset1);
    printf("\n");
    assert(vecset1.data[0].x == 0 && vecset1.data[0].y == 0);
    assert(vecset1.data[1].x == 5 && vecset1.data[1].y == 5);
    assert(vecset1.data[2].x == 10 && vecset1.data[2].y == 10);
    vecset_destroy(&vecset1);

    // Test vecset_push && vecset_top && vecset_second
    vec = vec_create(1, 1);
    vec1 = vec_create(2, 2);
    vecset_create(&vecset1);
    vecset_push(&vecset1, vec);
    assert(vecset1.data[0].x == 1 && vecset1.data[0].y == 1);
    vecset_push(&vecset1, vec1);
    assert(vecset1.data[0].x == 1 && vecset1.data[0].y == 1);
    assert(vecset1.data[1].x == 2 && vecset1.data[1].y == 2);
    const struct vec *top = vecset_top(&vecset1);
    assert(top->x == 2 && top->y == 2);
    const struct vec *second = vecset_second(&vecset1);
    assert(second->x == 1 && second->y == 1);

    // Test vecset_pop
    vecset_pop(&vecset1);
    vecset_push(&vecset1, vec);
    assert(vecset1.data[1].x == 1 && vecset1.data[1].y == 1);

    vecset_destroy(&vecset1);

    // Finished
    printf("Assertions passed");
}

int main() {
    tests();
    return 0;
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
