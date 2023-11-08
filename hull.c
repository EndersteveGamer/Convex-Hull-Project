// Théo Pariney

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stddef.h>

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

const int GROW_THRESHOLD = 2;

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
void vec_print(const struct vec *self) {
    printf("x: %f, y: %f", self->x, self->y);
}

/**
 * Returns the dot of two vectors
 * @param v1 The first vector
 * @param v2 The second vector
 * @return The dot of the two vectors
 */
double dot(const struct vec *v1, const struct vec *v2) {
    return v1->x * v2->x + v1->y * v2->y;
}

/**
 * Returns the cross product between three vectors
 * @param v1 The first vector
 * @param v2 The second vector
 * @param v3 The third vector
 * @return The cross product of the vectors
 */
double cross(const struct vec *v1, const struct vec *v2, const struct vec *v3) {
    return (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
}

/**
 * Checks if three vectors are in a left turn (cross product strictly negative)
 * @param v1 The first vector
 * @param v2 The second vector
 * @param v3 The third vector
 * @return true if the vectors are in a left turn, false otherwise
 */
bool is_left_turn(const struct vec *v1, const struct vec *v2, const struct vec *v3) {
    return cross(v1, v2, v3) > 0;
}

/**
 * Checks if two vectors are equal
 * @param v1 The first vector
 * @param v2 The second vector
 * @return true if the two vectors are equal, false otherwise
 */
bool vec_equals(const struct vec *v1, const struct vec *v2) {
    return v1->x == v2->x && v1->y == v2->y;
}

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
 * Prints a vecset in the same format than the hull-generator produces them
 * @param self The vecset to print
 */
void vecset_print(struct vecset *self) {
    printf("%zu\n", self->size);
    for (size_t i = 0; i < self->size; i++) {
        printf("%f %f\n", self->data[i].x, self->data[i].y);
    }
}

/**
 * Allocated more memory and increases the capacity of a vecset
 * @param self The vecset to grow
 */
void vecset_grow(struct vecset *self) {
    self->capacity *= GROW_THRESHOLD;
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

/**
 * Copies a vecset
 * @param self The vecset to copy
 * @return A pointer to the copied vecset
 */
struct vecset *vecset_copy(const struct vecset *self) {
    struct vecset *new = malloc(sizeof(struct vecset));
    vecset_create(new);
    for (size_t i = 0; i < self->size; i++) {
        struct vec newVec;
        newVec.x = self->data[i].x;
        newVec.y = self->data[i].y;
        vecset_add(new, newVec);
    }
    return new;
}

/**
 * A typedef corresponding to functions taking two vectors and comparing them, returning an int equal to 0 if they are
 * equal, inferior to 0 if the first vector is smaller than the second one, and superior to 0 otherwise
 */
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
    size_t index1 = 0;
    size_t index2 = 0;
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
    if (size <= 1) return;
    if (size == 2) {
        if (func(&array[0], &array[1], ctx) <= 0) return;
        struct vec temp = array[0];
        array[0] = array[1];
        array[1] = temp;
    }
    size_t size1 = size / 2;
    size_t size2 = size - size1;
    struct vec *a1 = calloc(size1, sizeof(struct vec));
    struct vec *a2 = calloc(size2, sizeof(struct vec));
    vec_array_split(a1, a2, array, size);
    vec_array_merge_sort(a1, size1, func, ctx);
    vec_array_merge_sort(a2, size2, func, ctx);
    vec_array_merge(a1, size1, a2, size2, array, func, ctx);
}

void vec_array_bubble_sort(struct vec *data, size_t n, comp_func_t func, const void *ctx) {
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = n - 1; j > i; --j) {
            if (func(&data[j], &data[j - 1], ctx) < 0) {
                struct vec temp = data[j];
                data[j] = data[j - 1];
                data[j - 1] = temp;
            }
        }
    }
}

/**
 * Sorts a vecset using a merge sort, according to a comparison function
 * @param self The vecset to sort
 * @param func The comparison function
 * @param ctx The context for the comparison function
 */
void vecset_sort(struct vecset *self, comp_func_t func, const void *ctx) {
    if (self->size <= 20) vec_array_bubble_sort(self->data, self->size, func, ctx);
    else vec_array_merge_sort(self->data, self->size, func, ctx);
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
 * <p>
 * Compares two points based on their distance to origin
 * </p>
 * <p>
 * Used only in unary tests
 * </p>
 * @param p1 The first point
 * @param p2 The second point
 * @param ctx The context for the function (unused)
 * @return A strictly positive number if p1 is further to the origin than p2, strictly negative if p2 is further, and
 * equal to 0 if both points have the same distance to the origin
 */
int comp_distance_to_origin(const struct vec *p1, const struct vec *p2, const void *ctx) {
    return (int)(p1->x + p1->y - p2->x - p2->y);
}

/**
 * Returns the leftmost point in a vecset (the point with the lowest x)
 * @param self The vecset to find the leftmost point in
 * @return The leftmost point in the vecset
 */
struct vec *vecset_leftmost_point(const struct vecset *self) {
    assert(self->size != 0);
    if (self->size == 1) return &self->data[0];
    struct vec *leftmost = self->data;
    struct vec *curr = self->data + 1;
    for (size_t i = 0; i < self->size - 1; i++) {
        if (curr->x < leftmost->x) leftmost = curr;
        curr++;
    }
    return leftmost;
}

/**
 * Given a vecset and a vector, finds the index of the vector in the vecset
 * @param self The vecset to search the vector in
 * @param vec The vector to find
 * @return The index of the vector if it is in the vecset, -1 otherwise
 */
size_t vecset_vec_index(const struct vecset *self, const struct vec *vec) {
    size_t index = 0;
    while (!vec_equals(&self->data[index], vec) && index < self->size) index++;
    return index < self->size ? index : -1;
}

/**
 * Unary tests on all the little functions (tests on convex hull algorithms are made manually)
 */
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

    // Test vecset_leftmost_point
    vecset_create(&vecset1);
    vec1 = vec_create(0, 0);
    vecset_add(&vecset1, vec3);
    vecset_add(&vecset1, vec2);
    vecset_add(&vecset1, vec1);
    struct vec *leftmost = vecset_leftmost_point(&vecset1);
    assert(leftmost->x == vec1.x && leftmost->y == vec1.y);
    vecset_destroy(&vecset1);

    // Finished
    printf("Assertions passed");
}

/**
 * A Jarvis March algorithm for finding the convex hull of a cloud of points
 * @param in The cloud of points to search the convex hull in
 * @param out A vecset containing all the points that are part of the convex hull
 */
void jarvis_march(const struct vecset *in, struct vecset *out) {
    if (in->size == 0) return;
    struct vec *F = vecset_leftmost_point(in);
    size_t FIndex = vecset_vec_index(in, F);
    struct vec *C = F;
    do {
        vecset_add(out, *C);
        FIndex = (FIndex + 1) % in->size;
        struct vec *N = &in->data[FIndex];
        for (size_t i = 0; i < in->size; i++) {
            if (i == FIndex) continue;
            if (is_left_turn(C, &in->data[i], N) || vec_equals(N, C)) N = &in->data[i];
        }
        C = N;
    } while (!vec_equals(F, C));
}

/**
 * Compared vectors based on their height. If both vectors have the same height, the one with the lowest absciss is
 * selected
 * @param v1 The first vector
 * @param v2 The second vector
 * @param ctx The context (unused)
 * @return 0 if both vectors are at the same location, 1 if the first vector is higher that the second, -1 if the second
 * vector is higher than the first
 */
int comp_vec_height(const struct vec *v1, const struct vec *v2, const void *ctx) {
    if (v1->y == v2->y) return v1->x == v2->x ? 0 : v1->x < v2->x ? 1 : -1;
    return v1->y < v2->y ? 1 : -1;
}

/**
 * Calculates the distance between two vectors
 * @param v1 The first vector
 * @param v2 The second vector
 * @return The distance between the two vectors
 */
double vec_distance(const struct vec *v1, const struct vec *v2) {
    return sqrt((v2->x - v1->x) * (v2->x - v1->x) + (v2->y - v1->y) * (v2->y - v1->y));
}

/**
 * <p>
 *      Compares vectors based on their angle with the vector in the context
 * </p>
 *  <p>
 *      If both angles are equals, the vectors are compared based on their distance with the vector in the context
 *  </p>
 * @param v1 The first vector
 * @param v2 The second vector
 * @param ctx A pointer to a vector
 * @return -1 if v1 has a lower angle, 1 if v2 has a lower angle, 1 if both vectors have the same angle and v1 is closer
 * to ctx than v2, else -1
 */
int comp_vec_angle(const struct vec *v1, const struct vec *v2, const void *ctx) {
    struct vec *B = (struct vec *) ctx;
    double v1B = atan2(B->y - v1->y, B->x - v1->x);
    double v2B = atan2(B->y - v2->y, B->x - v2->x);
    double difference = v1B - v2B;
    if (difference == 0) return vec_distance(B, v1) > vec_distance(B, v2) ? -1 : 1;
    return difference == 0 ? 0 : (difference > 0 ? 1 : -1);
}

void graham_scan(const struct vecset *in, struct vecset *out) {
    struct vecset *input = vecset_copy(in);
    const struct vec *B = vecset_min(in, &comp_vec_height, NULL);
    vecset_push(out, *B);
    vecset_sort(input, &comp_vec_angle, B);
    struct vec F = input->data[1];
    vecset_push(out, F);
    vecset_push(out, input->data[2]);
    for (size_t i = 3; i < input->size; i++) {
        if (vec_equals(&input->data[i], B)) continue;
        if (vec_equals(&input->data[i], &F)) continue;
        const struct vec *T = vecset_top(out);
        const struct vec *S = vecset_second(out);
        while (out->size >= 2 && !is_left_turn(S, T, &input->data[i])) {
            vecset_pop(out);
            T = vecset_top(out);
            S = vecset_second(out);
        }
        vecset_push(out, input->data[i]);
    }
    vecset_destroy(input);
    free(input);
}

/**
 * Returns the distance of the point M from the line XY using the equation to find the height of a triangle (with XY as
 * the base of the triangle)
 * @param X The first point of the line
 * @param Y The second point of the line
 * @param M The point to calculate the distance from the line
 * @return The distance of the point M from the line XY
 */
double vec_distance_to_line(const struct vec *X, const struct vec *Y, const struct vec *M) {
    double a = vec_distance(X, M);
    double b = vec_distance(X, Y);
    double c = vec_distance(Y, M);
    return (0.5 / b) * sqrt(a + b + c) * sqrt(-a + b + c) * sqrt(a - b + c) * sqrt(a + b - c);
}

/**
 * The second part of the quickhull algorithm
 * @param S A vecset
 * @param X The X vector of the quickhull algorithm
 * @param Y The Y vector of the quickhull algorithm
 * @return A vecset corresponding to the result of findHull for vectors on the left of the line XY and the result for
 * vectors on the right of the line XY
 */
struct vecset *findhull(const struct vecset *S, const struct vec *X, const struct vec *Y) {
    struct vecset *out = malloc(sizeof(struct vecset));
    vecset_create(out);
    if (S->size == 0) return out;
    struct vec *M = &S->data[0];
    double M_distance = vec_distance_to_line(X, Y, M);
    for (size_t i = 1; i < S->size; i++) {
        double i_distance = vec_distance_to_line(X, Y, &S->data[i]);
        if (i_distance > M_distance) {
            M = &S->data[i];
            M_distance = i_distance;
        }
    }
    struct vecset *s1 = malloc(sizeof(struct vecset));
    vecset_create(s1);
    struct vecset *s2 = malloc(sizeof(struct vecset));
    vecset_create(s2);
    for (size_t i = 0; i < S->size; i++) {
        if (vec_equals(M, &S->data[i])) continue;
        if (is_left_turn(X, M, &S->data[i])) vecset_add(s1, S->data[i]);
        if (is_left_turn(M, Y, &S->data[i])) vecset_add(s2, S->data[i]);
    }
    struct vecset *r1 = findhull(s1, X, M);
    vecset_destroy(s1);
    free(s1);
    struct vecset *r2 = findhull(s2, M, Y);
    vecset_destroy(s2);
    free(s2);
    for (size_t i = 0; i < r1->size; i++) vecset_add(out, r1->data[i]);
    vecset_add(out, *M);
    for (size_t i = 0; i < r2->size; i++) vecset_add(out, r2->data[i]);
    vecset_destroy(r1);
    free(r1);
    vecset_destroy(r2);
    free(r2);
    return out;
}

/**
 * Finds the convex hull of a vecset using a quickhull algorithm
 * @param in The vecset containing the vectors
 * @param out A vecset containing the vectors belonging to the convex hull
 */
void quickhull(const struct vecset *in, struct vecset *out) {
    struct vec *A = &in->data[0];
    struct vec *B = &in->data[0];
    // Find leftmost and rightmost point
    for (size_t i = 0; i < in->size; i++) {
        if (in->data[i].x < A->x) A = &in->data[i];
        if (in->data[i].x > B->x) B = &in->data[i];
    }
    struct vecset *s1 = malloc(sizeof(struct vecset));
    vecset_create(s1);
    struct vecset *s2 = malloc(sizeof(struct vecset));
    vecset_create(s2);
    for (size_t i = 0; i < in->size; i++) {
        if (vec_equals(&in->data[i], A)) continue;
        if (vec_equals(&in->data[i], B)) continue;
        if (is_left_turn(A, B, &in->data[i])) vecset_add(s1, in->data[i]);
        else vecset_add(s2, in->data[i]);
    }
    struct vecset *r1 = findhull(s1, A, B);
    vecset_destroy(s1);
    free(s1);
    struct vecset *r2 = findhull(s2, B, A);
    vecset_destroy(s2);
    free(s2);
    vecset_add(out, *A);
    for (size_t i = 0; i < r1->size; i++) vecset_add(out, r1->data[i]);
    vecset_add(out, *B);
    for (size_t i = 0; i < r2->size; i++) vecset_add(out, r2->data[i]);
    vecset_destroy(r1);
    free(r1);
    vecset_destroy(r2);
    free(r2);
}

/**
 * The main function executing the pilot for getting the file containing the cloud of points, parsing it, and printing
 * back the convex hull
 * @return 0 if everything executed correctly, or the error code otherwise
 */
int main() {
    // tests();
    // return 0;

    setbuf(stdout, NULL); // Avoid buffering in the output

    char buffer[BUFFSIZE];
    fgets(buffer, BUFFSIZE, stdin);

    size_t count = strtol(buffer, NULL, 10);

    struct vecset in;
    vecset_create(&in);
    struct vecset out;
    vecset_create(&out);

    for (size_t i = 0; i < count; ++i) {
        struct vec p;

        fgets(buffer, BUFFSIZE, stdin);

        char *endptr = buffer;
        p.x = strtod(endptr, &endptr);
        p.y = strtod(endptr, &endptr);

        vecset_add(&in, p);
    }

    jarvis_march(&in, &out);

    vecset_print(&out);

    vecset_destroy(&in);
    vecset_destroy(&out);
    return 0;
}
