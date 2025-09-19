#ifndef HASH_H
#define HASH_H

#define MAX_SIZE 100

typedef struct TypeInfo {

    size_t size;

    // Functions that define how to handle this data type
    int (*cmp_func)(const void*, const void*);
    size_t (*hash_func)(const void*);
    void (*free_func)(void*);
    void (*print_func)(const void*);

} TypeInfo;

// Linked list to store multiple values per key
typedef struct ValueNode {

    void* value;
    TypeInfo valueInfo;
    struct ValueNode* next;

} ValueNode;

// Linked list to store the keys
typedef struct KeyNode {

    void* key;
    TypeInfo keyInfo;
    ValueNode* val_head;
    struct KeyNode* next;

} KeyNode;

// Initalize the type information of a data type
void initType(TypeInfo* typeInfo, size_t size, int (*cmp_func)(const void*, const void*), size_t (*hash_func)(const void*),
              void (*free_func)(void*), void (*print_func)(const void*));

void hashtableAdd(KeyNode* table[], void* key_in, TypeInfo key_info, void* value_in, TypeInfo val_info);
void hashtableDelete(KeyNode* table[], void* key, TypeInfo key_info);
ValueNode* hashtableSearch(KeyNode* table[], void* key, TypeInfo key_info);
void hashtablePrint(KeyNode* table[]);
void hashtableFree(KeyNode* table[]);

#endif