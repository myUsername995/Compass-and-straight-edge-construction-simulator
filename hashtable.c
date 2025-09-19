#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include "hashtable.h"

void initType(TypeInfo* typeInfo, size_t size, int (*cmp_func)(const void*, const void*), size_t (*hash_func)(const void*),
              void (*free_func)(void*), void (*print_func)(const void*)){

    typeInfo->cmp_func = cmp_func;
    typeInfo->hash_func = hash_func;
    typeInfo->free_func = free_func;
    typeInfo->print_func = print_func;
    typeInfo->size = size;
}

// Use void value so that you can specify data type
void hashtableAdd(KeyNode* table[], void* key_in, TypeInfo key_info, void* value_in, TypeInfo val_info){
    int key_hashed = key_info.hash_func(key_in) % MAX_SIZE;

    KeyNode* node = table[key_hashed];
    KeyNode* prev = NULL;
    KeyNode* cur_key = NULL;

    while (node != NULL && key_info.cmp_func(node->key, key_in) != 0) {
        prev = node;
        node = node->next;
    }

    if (node != NULL) {
        cur_key = node;
    }

    // Create and initialize new value node
    ValueNode* new_val = malloc(sizeof(ValueNode));
    new_val->value = malloc(val_info.size);
    memcpy(new_val->value, value_in, val_info.size);
    new_val->valueInfo = val_info;
    new_val->next = NULL;

    if (!cur_key) {
        cur_key = malloc(sizeof(KeyNode));
        cur_key->key = malloc(key_info.size);
        memcpy(cur_key->key, key_in, key_info.size);
        cur_key->keyInfo = key_info;
        cur_key->val_head = new_val;
        cur_key->next = NULL;

        if (prev) {
            prev->next = cur_key;
        } else {
            table[key_hashed] = cur_key;
        }
    } else {
        ValueNode* temp = cur_key->val_head;
        while (temp->next){
            temp = temp->next;
        }
        temp->next = new_val;
    }
}

void hashtableDelete(KeyNode* table[], void* key, TypeInfo key_info){
    int key_hashed = key_info.hash_func(key) % MAX_SIZE;

    KeyNode* head = table[key_hashed];
    KeyNode* prev = NULL;
    KeyNode* cur = NULL;

    // Find the key node
    while (head != NULL) {
        if (key_info.cmp_func(key, head->key) == 0) {
            cur = head;
            break;
        }
        prev = head;
        head = head->next;
    }

    if (!cur){
        return;
    } // Key not found

    // Free all value nodes
    ValueNode* val = cur->val_head;
    while (val != NULL) {
        ValueNode* temp = val;
        val = val->next;

        if (temp->valueInfo.free_func) {
            temp->valueInfo.free_func(temp->value);
        } else {
            free(temp->value);
        }

        free(temp);
    }

    // Remove the key node from the list
    if (!prev) {
        table[key_hashed] = cur->next;
    } else {
        prev->next = cur->next;
    }

    if (key_info.free_func) {
        key_info.free_func(cur->key);
    } else {
        free(cur->key);
    }

    free(cur);
}

ValueNode* hashtableSearch(KeyNode* table[], void* key, TypeInfo key_info){
    int key_hashed = key_info.hash_func(key) % MAX_SIZE;

    KeyNode* node = table[key_hashed];

    while (node != NULL) {
        if (key_info.cmp_func(key, node->key) == 0){
            return node->val_head;
        }
        node = node->next;
    }

    //printf("Key not found!\n");
    return NULL; // Key not found
}

void hashtablePrint(KeyNode* table[]){
    for (size_t i = 0; i < MAX_SIZE; i++) {
        KeyNode* key_node = table[i];
        if (!key_node) continue;

        printf("Bucket %zu:\n", i);
        while (key_node) {
            printf("  Key: ");
            if (key_node->keyInfo.print_func) {
                key_node->keyInfo.print_func(key_node->key);
            } else {
                printf("[unknown type]");
            }
            printf(" -> Values: ");

            ValueNode* val_node = key_node->val_head;
            while (val_node) {
                if (val_node->valueInfo.print_func) {
                    val_node->valueInfo.print_func(val_node->value);
                } else {
                    printf("[unknown type]");
                }
                printf(" -> ");
                val_node = val_node->next;
            }
            printf("NULL\n");

            key_node = key_node->next;
        }
    }
}

// Free the values inside the hashtable, but dont free the original array (assume user is just reusing the hashtable)
void hashtableFree(KeyNode* table[]){
    for (int i = 0; i < MAX_SIZE; i++){
        // Free the keys
        // Empty slot
        if (!table[i]){
            continue;
        }

        KeyNode* keys = table[i];

        while (keys){
            // Free the values
            ValueNode* values = keys->val_head;

            while (values){
                // Free the value
                values->valueInfo.free_func(values->value);

                // Go onto the next one
                ValueNode* valTemp = values;
                values = values->next;

                free(valTemp);
            }

            // Free the key
            keys->keyInfo.free_func(keys->key);

            // Go onto the next one
            KeyNode* keyTemp = keys;
            keys = keys->next;

            free(keyTemp);
        }
    }

    // Set the contents of the table to NULL
    for (int i = 0; i < MAX_SIZE; i++){
        table[i] = NULL;
    }
}