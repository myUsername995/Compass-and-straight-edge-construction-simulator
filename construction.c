#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include "hashtable.h"

#define WINDOW_WIDTH 1000
#define WINDOW_HEIGHT 800

#define SAVE_FOLDER "saveFiles"
#define MAX_FILES 100

bool debug = false;

#define PI 3.14159

typedef enum State {
    DRAW_LINE,
    USE_COMPASS,
    NUM_ELEMENTS
} State;

typedef enum Dir {CLOCKWISE, ANTI_CLOCKWISE} Dir;

typedef enum Properties {USING, RADIUS, LOCKED, MOVING, ZOOM, HELP_LINE, FILL, CLR_TYPE} Properties;

typedef struct Array {
    void* arr;            // Just store as a void pointer as we will cast it anyway
    int track, size;
} Array;

typedef struct Line {

    SDL_FPoint p1;
    SDL_FPoint p2;
    bool helper_line;

} Line;

typedef struct AdjLines {

    SDL_FPoint p1;
    SDL_FPoint* adjPoints;
    int track, size;

} AdjLines;

typedef struct Circle {

    SDL_FPoint mid_p;
    float radius;
    double start_angle, end_angle;
    Dir rotate;
    bool helper_line;

} Circle;

// Structure that holds an array of intersection angles.
typedef struct {
    double angles[10];    // Maximum possible intersections is 8, but account for start and end angle.
    int count;
} AngleList;

typedef struct {
    double absAngle; // absolute angle (in degrees)
    double relAngle; // relative angle = absAngle - circle.start_angle (normalized to [0,360))
} AnglePair;

typedef struct Tex {

    SDL_Texture* tex;
    SDL_FRect pos;
    Properties props;

} Tex;

// Keep track of which regions should be filled
typedef struct fillRegions {

    SDL_Color color;
    SDL_FPoint clicked;             // Stores the absolute position of the click (not the screen position)

} fillRegions;

// idk goofy ahh struct name aaaa
typedef struct sortByDistance {

    SDL_FPoint intersect;
    double distance;

} sortByDistance;

typedef struct angleIndex {

    double angle;
    int index;

} angleIndex;

typedef struct linePoints {

    SDL_FPoint* points;
    int track, size;

} linePoints;

// Used to store all the data about the current construction
typedef struct SaveData {
    char name[256];             // The name of the file

    // The data stored in the file
	Array points;
	Array screen_points;
	Array lines;
	Array screen_lines;
	Array circles;
	Array screen_circles;
    Array fillRegions;

	double zoom;
	SDL_FPoint top_left;

} SaveData;

bool isBetween(float, float, float, Dir);
bool line_line_intersect(Line, Line, SDL_FPoint*);

// Create the functions for the hash table
int cmp_point(const void* a, const void* b){
    const SDL_FPoint* pa = (const SDL_FPoint*)a;
    const SDL_FPoint* pb = (const SDL_FPoint*)b;

    if (fabs(pa->x - pb->x) < 1e-2 && fabs(pa->y - pb->y) < 1e-2){
        return 0;
    }
    else {
        return -1;
    }
}

// Hashes based on where the point is
size_t hash_point(const void* a){
    const SDL_FPoint* pa = (const SDL_FPoint*)a;

    return pa->x + pa->y;
}

void free_point(void* a){
    SDL_FPoint* pa = (SDL_FPoint*)a;

    free(pa);
}

void print_point(const void* a){
    const SDL_FPoint* pa = (const SDL_FPoint*)a;

    printf("%f, %f\n", pa->x, pa->y);
}

int cmp_adj(const void* a, const void* b){
    return 0;
}

size_t hash_adj(const void* a){
    return 0;
}

void free_adj(void* a){
    AdjLines* pa = (AdjLines*)a;

    free(pa->adjPoints);
    free(pa);
}

void print_adj(const void* a){
    const AdjLines* pa = (const AdjLines*)a;

    printf("Point: %f, %f ", pa->p1.x, pa->p1.y);
    printf("Adjacent points: ");
    for (int i = 0; i < pa->track; i++){
        printf("%f, %f ", pa->adjPoints[i].x, pa->adjPoints[i].y);
    }
}

// Sort by the X values
int qsortCompPoint(const void* a, const void* b){
    const SDL_FPoint* pa = (const SDL_FPoint*)a;
    const SDL_FPoint* pb = (const SDL_FPoint*)b;

    return pa->x - pb->x;
}

int qsortCompareDistances(const void* a, const void* b){
    const sortByDistance* pa = (const sortByDistance*)a;
    const sortByDistance* pb = (const sortByDistance*)b;

    return pa->distance - pb->distance;
}

int qsortCompareAngle(const void* a, const void* b){
    const angleIndex* pa = (const angleIndex*)a;
    const angleIndex* pb = (const angleIndex*)b;

    return pa->angle - pb->angle;
}

int comparePoints(const SDL_FPoint* a, const SDL_FPoint* b) {
    if (a->x != b->x)
        return (a->x < b->x) ? -1 : 1;
    if (a->y != b->y)
        return (a->y < b->y) ? -1 : 1;
    return 0;
}

int qsortCompareLines(const void* a, const void* b) {
    const Line* la = (const Line*)a;
    const Line* lb = (const Line*)b;

    // Normalize endpoints for line A
    const SDL_FPoint* a1 = &la->p1;
    const SDL_FPoint* a2 = &la->p2;
    if (comparePoints(a2, a1) < 0) {
        const SDL_FPoint* temp = a1;
        a1 = a2;
        a2 = temp;
    }

    // Normalize endpoints for line B
    const SDL_FPoint* b1 = &lb->p1;
    const SDL_FPoint* b2 = &lb->p2;
    if (comparePoints(b2, b1) < 0) {
        const SDL_FPoint* temp = b1;
        b1 = b2;
        b2 = temp;
    }

    // Compare normalized endpoints
    int cmp1 = comparePoints(a1, b1);
    if (cmp1 != 0) return cmp1;

    return comparePoints(a2, b2);
}


void freeData(SaveData* dataArray[], int count){
    for (int i = 0; i < count; i++) {
        free(dataArray[i]->points.arr);
        free(dataArray[i]->screen_points.arr);
        free(dataArray[i]->lines.arr);
        free(dataArray[i]->screen_lines.arr);
        free(dataArray[i]->circles.arr);
        free(dataArray[i]->screen_circles.arr);
        free(dataArray[i]->fillRegions.arr);
        free(dataArray[i]);
    }
}

// Helper to write an Array
void write_array(FILE* fp, Array* array, int element_size) {
    fwrite(&array->track, sizeof(int), 1, fp);
    fwrite(&array->size, sizeof(int), 1, fp);
    fwrite(array->arr, element_size, array->track, fp);
}

// Helper to read an Array
void read_array(FILE* fp, Array* array, int element_size) {
    fread(&array->track, sizeof(int), 1, fp);
    fread(&array->size, sizeof(int), 1, fp);
    array->arr = malloc(element_size * array->track);
    fread(array->arr, element_size, array->track, fp);
}

char* findFile(const char inputName[256]) {
    WIN32_FIND_DATA findData;
    HANDLE hFind;

    char searchPath[MAX_PATH];
    snprintf(searchPath, MAX_PATH, "saveFiles\\*");

    hFind = FindFirstFile(searchPath, &findData);
    if (hFind == INVALID_HANDLE_VALUE) {
        return NULL;
    }

    char* result = NULL;
    do {
        if (!(findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)){
            char text[300];
            snprintf(text, sizeof(text), "%s.sav", inputName);
            if (strcmp(findData.cFileName, text) == 0) {
                // Allocate and return full path
                result = malloc(MAX_PATH);
                snprintf(result, MAX_PATH, "saveFiles\\%s", text);
                break;
            }
        }
    } while (FindNextFile(hFind, &findData) != 0);

    FindClose(hFind);
    return result;  // NULL if not found
}

// Delete a file with the given name (without extension)
int deleteFile(const char inputName[256]) {
    char fullPath[MAX_PATH];
    snprintf(fullPath, sizeof(fullPath), "%s\\%s.sav", SAVE_FOLDER, inputName);

    if (DeleteFile(fullPath)) {
        return 1;
    } else {
        return 0;
    }
}

void saveAFile(SaveData* fileData){

    // Create a file to save to
    char path[256];
    snprintf(path, sizeof(path), "%s/%s.sav", SAVE_FOLDER, fileData->name);
    FILE* fp = fopen(path, "wb");

    // Write the name
    fwrite(&fileData->name, sizeof(fileData->name), 1, fp);

    // Write the arrays
    write_array(fp, &fileData->points, sizeof(SDL_FPoint));
    write_array(fp, &fileData->screen_points, sizeof(SDL_FPoint));
    write_array(fp, &fileData->lines, sizeof(Line));
    write_array(fp, &fileData->screen_lines, sizeof(Line));
    write_array(fp, &fileData->circles, sizeof(Circle));
    write_array(fp, &fileData->screen_circles, sizeof(Circle));
    write_array(fp, &fileData->fillRegions, sizeof(fillRegions));

    // All the other data
    fwrite(&fileData->zoom, sizeof(fileData->zoom), 1, fp);
    fwrite(&fileData->top_left, sizeof(fileData->top_left), 1 , fp);

    // Free the allocated memory
    free(fileData->points.arr);
    free(fileData->screen_points.arr);
    free(fileData->lines.arr);
    free(fileData->screen_lines.arr);
    free(fileData->circles.arr);
    free(fileData->screen_circles.arr);
    free(fileData->fillRegions.arr);

    fclose(fp);
}

// Find the file with input name, then overwrite that file
void saveToFile(SaveData* fileData, char inputName[256]) {
    char* path = findFile(inputName);
    if (!path){
        printf("Couldn't find the file!");
        return;
    }

    FILE* fp = fopen(path, "wb");
    if (!fp) {
        perror("Failed to open file for writing");
        free(path);
        return;
    }

    fwrite(fileData->name, sizeof(fileData->name), 1, fp);

    write_array(fp, &fileData->points, sizeof(SDL_FPoint));
    write_array(fp, &fileData->screen_points, sizeof(SDL_FPoint));
    write_array(fp, &fileData->lines, sizeof(Line));
    write_array(fp, &fileData->screen_lines, sizeof(Line));
    write_array(fp, &fileData->circles, sizeof(Circle));
    write_array(fp, &fileData->screen_circles, sizeof(Circle));
    write_array(fp, &fileData->fillRegions, sizeof(fillRegions));

    fwrite(&fileData->zoom, sizeof(double), 1, fp);
    fwrite(&fileData->top_left, sizeof(SDL_FPoint), 1, fp);

    // Free the allocated memory
    free(fileData->points.arr);
    free(fileData->screen_points.arr);
    free(fileData->lines.arr);
    free(fileData->screen_lines.arr);
    free(fileData->circles.arr);
    free(fileData->screen_circles.arr);
    free(fileData->fillRegions.arr);

    fclose(fp);
    free(path);
}

// Load all SaveData files into an array
int loadAllSaveFiles(SaveData* saves[], int maxFiles) {
    WIN32_FIND_DATA findData;
    HANDLE hFind;

    char searchPath[MAX_PATH];
    snprintf(searchPath, MAX_PATH, "%s\\*", SAVE_FOLDER);

    hFind = FindFirstFile(searchPath, &findData);
    if (hFind == INVALID_HANDLE_VALUE) {
        perror("Failed to open save folder");
        return 0;
    }

    int count = 0;
    do {
        if (!(findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            char fullPath[MAX_PATH];
            snprintf(fullPath, MAX_PATH, "%s\\%s", SAVE_FOLDER, findData.cFileName);

            FILE* fp = fopen(fullPath, "rb");
            if (!fp) continue;

            SaveData* data = malloc(sizeof(SaveData));
            fread(data->name, sizeof(data->name), 1, fp);

            read_array(fp, &data->points, sizeof(SDL_FPoint));
            read_array(fp, &data->screen_points, sizeof(SDL_FPoint));
            read_array(fp, &data->lines, sizeof(Line));
            read_array(fp, &data->screen_lines, sizeof(Line));
            read_array(fp, &data->circles, sizeof(Circle));
            read_array(fp, &data->screen_circles, sizeof(Circle));
            read_array(fp, &data->fillRegions, sizeof(fillRegions));

            fread(&data->zoom, sizeof(double), 1, fp);
            fread(&data->top_left, sizeof(SDL_FPoint), 1, fp);

            fclose(fp);
            saves[count++] = data;

            if (count >= maxFiles) break;
        }
    } while (FindNextFile(hFind, &findData) != 0);

    FindClose(hFind);
    return count;
}

// Take an x,y position as input for position, and add width and height of the texture to it
// Also create a texture from just a string
SDL_Texture* renderTexture(SDL_Renderer* renderer, TTF_Font* font, char* str, SDL_FRect* pos, SDL_Color color){
    SDL_Surface* surface = TTF_RenderText_Solid(font, str, strlen(str), color);

    pos->w = surface->w;
    pos->h = surface->h;

    SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);

    SDL_DestroySurface(surface);

    return texture;
}

SDL_FPoint world_to_screen(SDL_FPoint p, double zoom, SDL_FPoint top_left){
    return (SDL_FPoint){(p.x - top_left.x) / zoom, (p.y - top_left.y) / zoom};
}

SDL_FPoint screen_to_world(SDL_FPoint p, double zoom, SDL_FPoint top_left){
    return (SDL_FPoint){p.x * zoom + top_left.x, p.y * zoom + top_left.y};
}

void drawPoint(SDL_Renderer* renderer, SDL_FPoint point, int radius){

    // Use scanline rendering for a circle
    for (int i = -radius; i < radius; i++){
        int offset = sqrt(radius * radius - i * i);

        SDL_RenderLine(renderer, point.x + offset, point.y - i, point.x - offset, point.y - i);
    }
}

double normalizeAngle(double angle) {
    return fmod(angle + 360.0, 360.0);
}

int comparePair(const void* a, const void* b) {
    const AnglePair* pa = (const AnglePair*)a;
    const AnglePair* pb = (const AnglePair*)b;

    if (pa->relAngle < pb->relAngle){
        return -1;
    }
    else if (pa->relAngle > pb->relAngle){
        return 1;
    }
    else {
        return 0;
    }
}

AngleList getIntersectionAngles(Circle circle) {
    AngleList angleList;
    angleList.count = 0;

    double cx = circle.mid_p.x;
    double cy = circle.mid_p.y;
    double r = circle.radius;

    bool full_circle = false;

    // Handle a near full circle (which messes up the rendering)
    if (fabs(circle.start_angle - circle.end_angle) < 1e-4){
        full_circle = true;
        if (circle.rotate == CLOCKWISE){
            circle.start_angle += 1e-4;
        }
        else if (circle.rotate == ANTI_CLOCKWISE){
            circle.start_angle -= 1e-4;
        }
    }

    // Optionally, add the circle's start angle if its point is visible.
    double start_x = cx + r * cos(circle.start_angle * PI / 180.0);
    double start_y = cy + r * sin(circle.start_angle * PI / 180.0);

    if (start_x >= 0 && start_x <= WINDOW_WIDTH && start_y >= 0 && start_y <= WINDOW_HEIGHT){
        angleList.angles[angleList.count++] = circle.start_angle;
    }

    // Check vertical edges: x = 0 and x = WINDOW_WIDTH.
    double xs[2] = {0.0, (double)WINDOW_WIDTH};
    for (int i = 0; i < 2; i++) {
        double xPos = xs[i];
        double dx = xPos - cx;
        double discriminant = r * r - dx * dx;
        if (discriminant >= 0) {
            double root = sqrt(discriminant);
            double y1 = cy + root;
            double y2 = cy - root;
            
            if (y1 >= 0 && y1 <= WINDOW_HEIGHT) {
                double angle = normalizeAngle(atan2(y1 - cy, dx) * 180.0 / PI);
                if (isBetween(circle.start_angle, circle.end_angle, angle, circle.rotate) || full_circle){
                    angleList.angles[angleList.count++] = angle;
                }
            }
            if (y2 >= 0 && y2 <= WINDOW_HEIGHT) {
                double angle = normalizeAngle(atan2(y2 - cy, dx) * 180.0 / PI);
                if (isBetween(circle.start_angle, circle.end_angle, angle, circle.rotate) || full_circle) {
                    angleList.angles[angleList.count++] = angle;
                }
            }
        }
    }
    
    // Check horizontal edges: y = 0 and y = WINDOW_HEIGHT.
    double ys[2] = {0.0, (double)WINDOW_HEIGHT};
    for (int i = 0; i < 2; i++) {
        double yPos = ys[i];
        double dy = yPos - cy;
        double discriminant = r * r - dy * dy;
        if (discriminant >= 0) {
            double root = sqrt(discriminant);
            double x1 = cx + root;
            double x2 = cx - root;
            
            if (x1 >= 0 && x1 <= WINDOW_WIDTH) {
                double angle = normalizeAngle(atan2(yPos - cy, x1 - cx) * 180.0 / PI);
                if (isBetween(circle.start_angle, circle.end_angle, angle, circle.rotate) || full_circle) {
                    angleList.angles[angleList.count++] = angle;
                }
            }
            if (x2 >= 0 && x2 <= WINDOW_WIDTH) {
                double angle = normalizeAngle(atan2(yPos - cy, x2 - cx) * 180.0 / PI);
                if (isBetween(circle.start_angle, circle.end_angle, angle, circle.rotate) || full_circle) {
                    angleList.angles[angleList.count++] = angle;
                }
            }
        }
    }

    // Optionally, add the circle's end angle if visible.
    double end_x = cx + r * cos(circle.end_angle * PI / 180.0);
    double end_y = cy + r * sin(circle.end_angle * PI / 180.0);

    if (end_x >= 0 && end_x <= WINDOW_WIDTH && end_y >= 0 && end_y <= WINDOW_HEIGHT) {
        angleList.angles[angleList.count++] = circle.end_angle;
    }

    double rel_angle = circle.start_angle;
    AnglePair pairs[10];

    // Use some goofy ahh function to calculate clockwise and anti-clockwise distance
    for (int i = 0; i < angleList.count; i++){
        if (circle.rotate == ANTI_CLOCKWISE){
            pairs[i].relAngle = fmod(rel_angle - angleList.angles[i] + 360.0, 360.0);
        }
        else {
            pairs[i].relAngle = fmod(angleList.angles[i] - rel_angle + 360.0, 360.0);
        }
        pairs[i].absAngle = angleList.angles[i];
    }

    // Sort relative angles
    qsort(pairs, angleList.count, sizeof(AnglePair), comparePair);

    // Copy the angles into the angleList.angles array
    for (int i = 0; i < angleList.count; i++){
        angleList.angles[i] = pairs[i].absAngle;
    }

    return angleList;
}

void drawSegmentCircle(SDL_Renderer* renderer, Circle circle) {
    bool clockwise = circle.rotate == CLOCKWISE;
    double start = circle.start_angle;
    double end = circle.end_angle;
    
    // Compute the arc length based on rotation direction.
    // If start == end, assume a full circle (360 degrees).
    float arcLength = 0.0f;
    if (clockwise) {
        // For clockwise, if start is less than end, add 360.
        arcLength = (start - end);
        if (arcLength <= 0)
            arcLength = (start == end ? 360.0f : (start - end + 360.0f));
    } else {
        // For counterclockwise, if end is less than start, add 360.
        arcLength = (end - start);
        if (arcLength <= 0)
            arcLength = (start == end ? 360.0f : (end - start + 360.0f));
    }
    
    // Choose a step; here we use a degree per step.
    int steps;
    if (circle.helper_line){
        steps = arcLength * circle.radius / 200;
    }
    else {
        steps = arcLength * circle.radius / 10;
    }
    if (steps < 1) steps = 1;
    float stepSize = arcLength / steps;

    SDL_FPoint points[2];
    int track = 0;
    bool first_point = true;

    // Draw the arc
    for (int i = 0; i <= steps; i++) {
        float angle;
        if (clockwise) {
            angle = start - i * stepSize;
        } else {
            angle = start + i * stepSize;
        }
        // Normalize the angle to keep it within [0, 360)
        angle = fmod(angle + 360.0, 360.0f);

        double x = circle.mid_p.x + circle.radius * cos(angle * PI / 180.0);
        double y = circle.mid_p.y + circle.radius * sin(angle * PI / 180.0);

        points[track] = (SDL_FPoint){x, y};

        track++;
        track %= 2;

        if (!first_point){
            SDL_RenderLine(renderer, points[track].x, points[track].y, points[(track+1)%2].x, points[(track+1)%2].y);
        }
        else {
            first_point = false;
        }
    }
}

void drawClippedSegmentCircles(SDL_Renderer* renderer, Circle circle) {
    AngleList list = getIntersectionAngles(circle);

    for (int i = 0; i < list.count; i += 2){
        circle.end_angle = list.angles[i];
        circle.start_angle = list.angles[i+1];

        drawSegmentCircle(renderer, circle);
    }
}

void drawClippedLine(SDL_Renderer* renderer, Line line){
    // Calculate the intersections with the edges of the screen first

    SDL_FPoint intersects[2];            // Amount of intersections possible
    int track = 0;

    // The four edges of the screen
    Line sides[4] = {{.p1.x = 0, .p1.y = 0, .p2.x = 0, .p2.y = WINDOW_HEIGHT, .helper_line = false},
                     {.p1.x = 0, .p1.y = 0, .p2.x = WINDOW_WIDTH, .p2.y = 0, . helper_line = false},
                     {.p1.x = WINDOW_WIDTH, .p1.y = 0, .p2.x = WINDOW_WIDTH, .p2.y = WINDOW_HEIGHT},
                     {.p1.x = 0, .p1.y = WINDOW_HEIGHT, .p2.x = WINDOW_WIDTH, .p2.y = WINDOW_HEIGHT}};

    for (int i = 0; i < 4; i++){
        SDL_FPoint intersect;

        // Function to check bounded line intersections
        bool isIntersected = line_line_intersect(line, sides[i], &intersect);

        if (isIntersected && track < 2){
            intersects[track++] = intersect;
        }
    }

    if (track > 2){
        return;
    }

    // Draw the line with the intersection points
    if (track == 0){
        SDL_RenderLine(renderer, line.p1.x, line.p1.y, line.p2.x, line.p2.y);
    }
    // Check which point is inside the visible space, then draw between that point and the intersection
    else if (track == 1){
        if (line.p1.x > 0 && line.p1.x < WINDOW_WIDTH && line.p1.y > 0 && line.p1.y < WINDOW_HEIGHT){
            SDL_RenderLine(renderer, line.p1.x, line.p1.y, intersects[0].x, intersects[0].y);
        }
        else {
            SDL_RenderLine(renderer, line.p2.x, line.p2.y, intersects[0].x, intersects[0].y);
        }
    }
    // Draw between the two intersections as it is the only part visible
    else if (track == 2){
        SDL_RenderLine(renderer, intersects[0].x, intersects[0].y, intersects[1].x, intersects[1].y);
    }
}

void calculateSlope(float x1, float y1, float x2, float y2, double* slope, double* intercept){
    if (x2 - x1 == 0){
        *slope = INFINITY;

        return;
    }
    *slope = (y2 - y1) / (x2 - x1);
    *intercept = y2 - x2 * (*slope);
}

bool isBetween(float start, float end, float mid, Dir direction) {     
    end = (end - start) < 0.0f ? end - start + 360.0f : end - start;    
    mid = (mid - start) < 0.0f ? mid - start + 360.0f : mid - start;
    
    if (direction == CLOCKWISE){
        return (mid < end);
    }
    else {
        return (mid >= end);
    }
}

// Checks where two circles intersect, returns false when no intersection happens, true when it does
bool circle_circle_intersection(SDL_FPoint point1, double r1, SDL_FPoint point2, double r2, SDL_FPoint *p1, SDL_FPoint *p2){
    double x1 = point1.x;
    double y1 = point1.y;

    double x2 = point2.x;
    double y2 = point2.y;

    // Distance between centers
    double d = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Check if circles intersect
    if (d > r1 + r2 || d < fabs(r1 - r2)) {
        return false;  // No intersection
    }

    // Find midpoint along the line connecting centers
    double a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    double h = sqrt(r1 * r1 - a * a);

    // Coordinates of midpoint
    double x_mid = x1 + (a / d) * (x2 - x1);
    double y_mid = y1 + (a / d) * (y2 - y1);

    // Intersection points
    p1->x = x_mid + (h / d) * (y2 - y1);
    p1->y = y_mid - (h / d) * (x2 - x1);

    p2->x = x_mid - (h / d) * (y2 - y1);
    p2->y = y_mid + (h / d) * (x2 - x1);

    return true;       // Intersection happened
}

void circle_segment_intersection(Circle circle1, Circle circle2, SDL_FPoint intersections[2], int* track){
    SDL_FPoint p1, p2;

    bool intersected = circle_circle_intersection(circle1.mid_p, circle1.radius, circle2.mid_p, circle2.radius, &p1, &p2);

    if (fabs(circle1.start_angle - circle1.end_angle) < 1e-4){
        if (circle1.rotate == CLOCKWISE){
            circle1.start_angle += 1e-4;
        }
        else if (circle1.rotate == ANTI_CLOCKWISE){
            circle1.start_angle -= 1e-4;
        }
    }

    if (fabs(circle2.start_angle - circle2.end_angle) < 1e-4){
        if (circle2.rotate == CLOCKWISE){
            circle2.start_angle += 1e-4;
        }
        else if (circle2.rotate == ANTI_CLOCKWISE){
            circle2.start_angle -= 1e-4;
        }
    }

    // No intesection between the two circles
    if (!intersected){
        *track = 0;
        return;
    }

    *track = 0;

    // Angle calculation
    double dist_x = p1.x - circle1.mid_p.x;
    double dist_y = p1.y - circle1.mid_p.y;

    double angle = atan2(dist_y, dist_x) * 180.0 / PI;

    if (angle < 0){
        angle += 360.0;
    }

    // Check if the intersection point lies in the circle segment
    bool intersectsSegment = isBetween(circle1.start_angle, circle1.end_angle, angle, circle1.rotate);

    // Same logic but apply for the other circle
    dist_x = p1.x - circle2.mid_p.x;
    dist_y = p1.y - circle2.mid_p.y;

    angle = atan2(dist_y, dist_x) * 180.0 / PI;

    if (angle < 0){
        angle += 360.0;
    }

    // Check if the intersection point lies in the circle segment
    bool intersectsSegment2 = isBetween(circle2.start_angle, circle2.end_angle, angle, circle2.rotate);

    if (intersectsSegment && intersectsSegment2){
        intersections[*track] = p1;
        (*track)++;
    }

    // Angle calculation
    dist_x = p2.x - circle1.mid_p.x;
    dist_y = p2.y - circle1.mid_p.y;

    angle = atan2(dist_y, dist_x) * 180.0 / PI;

    if (angle < 0){
        angle += 360.0;
    }

    // Check if the intersection point lies in the circle segment
    intersectsSegment = isBetween(circle1.start_angle, circle1.end_angle, angle, circle1.rotate);

    // Same logic but apply for the other circle
    dist_x = p2.x - circle2.mid_p.x;
    dist_y = p2.y - circle2.mid_p.y;

    angle = atan2(dist_y, dist_x) * 180.0 / PI;

    if (angle < 0){
        angle += 360.0;
    }

    // Check if the intersection point lies in the circle segment
    intersectsSegment2 = isBetween(circle2.start_angle, circle2.end_angle, angle, circle2.rotate);

    if (intersectsSegment && intersectsSegment2){
        intersections[*track] = p2;
        (*track)++;
    }
}

void circle_line_intersect(Circle circle, Line line, SDL_FPoint arr[2], int* intersect_track){
    float x0 = circle.mid_p.x;
    float y0 = circle.mid_p.y;

    double slope, intercept, discriminant;

    SDL_FPoint* intersections = arr;

    calculateSlope(line.p1.x, line.p1.y, line.p2.x, line.p2.y, &slope, &intercept);

    if (slope == INFINITY){
        discriminant = circle.radius * circle.radius - (line.p1.x - x0) * (line.p1.x - x0);

        if (discriminant < 0){
            *intersect_track = 0;
        }
        else if (discriminant == 0){
            intersections[0] = (SDL_FPoint){line.p1.x, y0};

            *intersect_track = 1;
        }
        else {
            float y1 = y0 + sqrt(discriminant);
            float y2 = y0 - sqrt(discriminant);

            intersections[0] = (SDL_FPoint){line.p1.x, y1};
            intersections[1] = (SDL_FPoint){line.p1.x, y2};

            *intersect_track = 2;
        }
    }
    else {
        double a = 1.0 + slope * slope;
        double b = 2.0 * (slope * (intercept - y0) - x0);
        double c = x0 * x0 + (intercept - y0) * (intercept - y0) - circle.radius * circle.radius;

        discriminant = b * b - 4 * a * c;

        if (discriminant < 0){
            *intersect_track = 0;
        }
        else if (discriminant == 0){
            float x = (-b) / (2.0f * a);
            float y = slope * x + intercept;

            intersections[0] = (SDL_FPoint){x, y};

            *intersect_track = 1;
        }
        else {
            float x1 = (-b + sqrt(discriminant)) / (2 * a);
            float y1 = slope * x1 + intercept;

            float x2 = (-b - sqrt(discriminant)) / (2 * a);
            float y2 = slope * x2 + intercept;

            intersections[0] = (SDL_FPoint){x1, y1};
            intersections[1] = (SDL_FPoint){x2, y2};

            *intersect_track = 2;
        }
    }

    // Check if the segmented line actually intersects the point
    SDL_FPoint new_intersect[2];
    int new_track = 0;

    for (int i = 0; i < *intersect_track; i++){
        bool bool_x1 = line.p1.x <= intersections[i].x && line.p2.x >= intersections[i].x;
        bool bool_y1 = line.p1.y <= intersections[i].y && line.p2.y >= intersections[i].y;

        bool bool_x2 = line.p1.x >= intersections[i].x && line.p2.x <= intersections[i].x;
        bool bool_y2 = line.p1.y >= intersections[i].y && line.p2.y <= intersections[i].y;

        if ((bool_x1 || bool_x2) && (bool_y1 || bool_y2)){
            new_intersect[new_track++] = intersections[i];
        }
    }

    int new_new_track = 0;

    // Check if the line segment intersects the circle segment
    for (int i = 0; i < new_track; i++){
        double dist_x = new_intersect[i].x - circle.mid_p.x;
        double dist_y = new_intersect[i].y - circle.mid_p.y;

        double angle = atan2(dist_y, dist_x) * 180.0 / PI;

        if (angle < 0){
            angle += 360.0;
        }

        if (fabs(circle.start_angle - circle.end_angle) < 1e-4){
            if (circle.rotate == CLOCKWISE){
                circle.start_angle += 1e-4;
            }
            else if (circle.rotate == ANTI_CLOCKWISE){
                circle.start_angle -= 1e-4;
            }
        }

        bool contains = isBetween(circle.start_angle, circle.end_angle, angle, circle.rotate);

        if (contains){
            new_intersect[new_new_track++] = new_intersect[i];
        }
    }

    memcpy(arr, new_intersect, new_new_track * sizeof(SDL_FPoint));
    *intersect_track = new_new_track;
}

bool line_line_intersect(Line line1, Line line2, SDL_FPoint* intersection){
    // Compute the denominator D
    double D = (line1.p1.x - line1.p2.x) * (line2.p1.y - line2.p2.y)
             - (line1.p1.y - line1.p2.y) * (line2.p1.x - line2.p2.x);

    // Parallel or colinear
    if (fabs(D) < 1e-9) {
        return false;
    }

    // Compute the parameters t and u
    double t = ((line1.p1.x - line2.p1.x) * (line2.p1.y - line2.p2.y)
              - (line1.p1.y - line2.p1.y) * (line2.p1.x - line2.p2.x)) 
              / D;

    double u = ((line1.p1.x - line2.p1.x) * (line1.p1.y - line1.p2.y)
              - (line1.p1.y - line2.p1.y) * (line1.p1.x - line1.p2.x)) 
              / D;

    // Check if intersection lies on both segments
    if (t < 0.0 || t > 1.0 || u < 0.0 || u > 1.0) {
        return false;
    }


    // Compute intersection point
    intersection->x = line1.p1.x + t * (line1.p2.x - line1.p1.x);
    intersection->y = line1.p1.y + t * (line1.p2.y - line1.p1.y);

    return true;
}

void drawCircle(SDL_Renderer* renderer, SDL_FPoint point, int radius){
    // Use the scanline but dont render filled circles

    Circle circle = {.mid_p = point, .radius = radius, .rotate = CLOCKWISE, .helper_line = false, .start_angle = 0, 
                     .end_angle = 0};

    AngleList angles = getIntersectionAngles(circle);

    // Needs to be clipped
    if (angles.count > 0){

        // For loop to draw the segment circles between the clipped points
        for (int i = 0; i < angles.count; i += 2){
            Circle new_circle = circle;

            new_circle.end_angle = angles.angles[i];
            new_circle.start_angle = angles.angles[i+1];

            drawSegmentCircle(renderer, new_circle);
        }

        // Dont continue with the function
        return;
    }

    SDL_FPoint left_side[2];
    SDL_FPoint right_side[2];

    int track = 0;
    for (int i = -radius; i <= radius; i++){
        int offset = sqrt(radius * radius - i * i);

        left_side[track] = (SDL_FPoint){point.x - offset, point.y + i};
        right_side[track] = (SDL_FPoint){point.x + offset, point.y + i};

        int next_index = (track + 1) % 2;
        int cur_index = track;

        track++;
        track %= 2;

        // Starting point (connect at the top)
        if (i == -radius){
            SDL_RenderLine(renderer, left_side[cur_index].x, left_side[cur_index].y, right_side[cur_index].x, 
                           right_side[cur_index].y);
            continue;
        }

        // Ending point (connect at the bottom)
        if (i == radius){
            SDL_RenderLine(renderer, left_side[next_index].x, left_side[next_index].y, right_side[next_index].x, 
                           right_side[next_index].y);
            continue;
        }

        // Connect adjacent points
        SDL_RenderLine(renderer, right_side[next_index].x, right_side[next_index].y, right_side[cur_index].x, 
                       right_side[cur_index].y);
        SDL_RenderLine(renderer, left_side[next_index].x, left_side[next_index].y, left_side[cur_index].x, 
                       left_side[cur_index].y);
    }
}

void resize(void** arr, int size, int elements_size){
    void* new_arr = realloc(*arr, size * elements_size);

    if (!new_arr){
        printf("Malloc fail\n");
        return;
    }

    *arr = new_arr;
}

SDL_FPoint nearest_point(SDL_FPoint* points, int track, float clicked_x, float clicked_y, float zoom){
    double min_point = INFINITY;
    SDL_FPoint nearest;

    for (int i = 0; i < track; i++){
        float dist_x = points[i].x - clicked_x;
        float dist_y = points[i].y - clicked_y;

        double dist = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);

        if (dist < min_point){
            min_point = dist;

            nearest = points[i];
        }
    }

    if (min_point < 30 * zoom){
        return nearest;
    }
    else {
        return (SDL_FPoint){clicked_x, clicked_y};
    }
}

void delete_arr(void *arr, int index, int count, int elem_size) {
    if (index < 0 || index >= count) {
        printf("Invalid index!\n");
        return;
    }
    
    char *char_arr = (char *)arr; // Treat array as char buffer

    // Use memmove (instead of memcpy) because src and dst overlap, which would cause memcpy to output undefined behaviour
    memmove(char_arr + index * elem_size, 
            char_arr + (index + 1) * elem_size, 
            (count - index - 1) * elem_size);
}

double crossProduct(double x1, double y1, double x2, double y2){
    return x1 * y2 - y1 * x2;
}

bool equalPoints(SDL_FPoint p1, SDL_FPoint p2){
    return fabs(p1.x - p2.x) < 1e-2 && fabs(p1.y - p2.y) < 1e-2;
}

// The order doesn't matter for the lines
bool equalLines(Line line1, Line line2){
    return equalPoints(line1.p1, line2.p1) && equalPoints(line1.p2, line2.p2) ||
           equalPoints(line1.p1, line2.p2) && equalPoints(line1.p2, line2.p1);
}

int getLineIndex(Line* lines, int line_track, Line line){
    for (int i = 0; i < line_track; i++){
        if (equalLines(lines[i], line)){
            return i;
        }
    }

    return 0;
}

void addPoint(KeyNode* hashtable[], SDL_FPoint point, SDL_FPoint adjPoint, TypeInfo pointInfo, TypeInfo adjInfo){
    // Add the first node of the line to the hashtable first
    ValueNode* hashSearch = hashtableSearch(hashtable, &point, pointInfo);
    // Check if the value exists in the hashtable
    if (!hashSearch){
        // Just add it to the hashtable in this case
        AdjLines adjList;
        adjList.p1 = point;
        adjList.adjPoints = malloc(100 * sizeof(SDL_FPoint));
        adjList.adjPoints[0] = adjPoint;
        adjList.track = 1;
        adjList.size = 100;

        hashtableAdd(hashtable, &adjList.p1, pointInfo, &adjList, adjInfo);
    }
    // Add it to an already existing value
    else {
        AdjLines adjList = *(AdjLines*)hashSearch->value;
        adjList.adjPoints[adjList.track++] = adjPoint;

        // Resize if necessary
        if (adjList.track >= adjList.size){
            resize((void**)&adjList.adjPoints, adjList.size, sizeof(SDL_FPoint));
        }

        // Re-alloc the lists
        SDL_FPoint* newPoints = malloc(adjList.size * sizeof(SDL_FPoint));
        memcpy(newPoints, adjList.adjPoints, adjList.size * sizeof(SDL_FPoint));

        // Delete it first
        hashtableDelete(hashtable, &adjList.p1, pointInfo);

        // Re-assing the list because it was freed in delete
        adjList.adjPoints = newPoints;
        
        // Then add it
        hashtableAdd(hashtable, &adjList.p1, pointInfo, &adjList, adjInfo);
    }
}

// Removes the element specified by adjPoint inside the adjacency list specified by point
void removePoint(KeyNode* hashtable[], SDL_FPoint point, SDL_FPoint adjPoint, TypeInfo pointInfo, TypeInfo adjInfo){
    ValueNode* hashSearch = hashtableSearch(hashtable, &point, pointInfo);

    // The point doesn't exist
    if (!hashSearch){
        return;
    }
    AdjLines adjList = *(AdjLines*)hashSearch->value;

    // Search for the value inside the adjacency list
    for (int i = 0; i < adjList.track; i++){
        // If we find the element, then we delete it from the array
        if (equalPoints(adjList.adjPoints[i], adjPoint)){
            delete_arr(adjList.adjPoints, i, adjList.track, sizeof(SDL_FPoint));
            adjList.track--;
        }
    }

    // Allocate a new array as it will be deleted during hashtableDelete operation
    SDL_FPoint* newPoints = malloc(adjList.size * sizeof(SDL_FPoint));
    memcpy(newPoints, adjList.adjPoints, adjList.size * sizeof(SDL_FPoint));

    // Replace the value
    hashtableDelete(hashtable, &point, pointInfo);

    adjList.adjPoints = newPoints;
    hashtableAdd(hashtable, &point, pointInfo, &adjList, adjInfo);
}

// Keeps the lines array input the same, and returns a new array

// Splits a line with n intersection points (not including endpoints) into n+1 segments
void splitLine(KeyNode* hashtable[], SDL_FPoint* intersects, int int_track, Line curLine, Line** newLines, 
               int* line_track, int* line_size, TypeInfo pointInfo, TypeInfo adjInfo){

    // Remove the main line
    for (int i = 0; i < *line_track; i++){
        if (equalPoints(curLine.p1, (*newLines)[i].p1) && equalPoints(curLine.p2, (*newLines)[i].p2) ||
            equalPoints(curLine.p2, (*newLines)[i].p1) && equalPoints(curLine.p1, (*newLines)[i].p2)){

            removePoint(hashtable, curLine.p1, curLine.p2, pointInfo, adjInfo);
            removePoint(hashtable, curLine.p2, curLine.p1, pointInfo, adjInfo);

            delete_arr(newLines, i, *line_track, sizeof(Line));
            (*line_track)--;
        }
    }

    // Compute all distances from the first point of the array (line[cur].p1)
    // and sort them

    sortByDistance* distances = malloc(int_track * sizeof(sortByDistance));
    distances[0].intersect = intersects[0];
    distances[0].distance = 0;

    // Find the distances from the start of every point and attach the point itself
    // to the array of structs
    for (int i = 1; i < int_track; i++){
        distances[i].intersect = intersects[i];

        double dist_x = intersects[0].x - intersects[i].x;
        double dist_y = intersects[0].y - intersects[i].y;

        distances[i].distance = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);
    }

    // Sort this tuff array
    qsort(distances, int_track, sizeof(sortByDistance), qsortCompareDistances);

    // Connect the lines in a tuff way
    for (int i = 0; i < int_track-1; i++){

        addPoint(hashtable, distances[i].intersect, distances[i+1].intersect, 
                    pointInfo, adjInfo);

        addPoint(hashtable, distances[i+1].intersect, distances[i].intersect, 
                    pointInfo, adjInfo);

        (*newLines)[(*line_track)++] = (Line){.p1 = distances[i].intersect, 
                                           .p2 = distances[i+1].intersect, 
                                           .helper_line = curLine.helper_line};

        if (*line_track >= *line_size){
            (*line_size) *= 2;
            resize((void**)newLines, *line_size, sizeof(Line));
        }
    }

    free(distances);
}

// Assume the newLines array is the same size, has the same track, and contains the same points as the lines array
// Split the lines so none of them are intersecting anymore
void splitAllLines(KeyNode* hashtable[], SDL_FPoint* points, int points_track, Line** newLines, int* new_track, int* new_size,
                   Line* lines, Line* screen_lines, int line_track, int line_size, TypeInfo pointInfo, TypeInfo adjInfo){

    // Get rid of the hashtable and make a new one
    hashtableFree(hashtable);

    // Get rid of the no intersect array and make a new one
    *new_track = 0;
    linePoints* lineP = malloc(line_track * sizeof(linePoints));

    // Initalise the array
    for (int i = 0; i < line_track; i++){
        lineP[i].points = malloc(100 * sizeof(SDL_FPoint));
        lineP[i].track = 0;
        lineP[i].size = 100;

        // Add the starting points
        lineP[i].points[lineP[i].track++] = lines[i].p1;
        lineP[i].points[lineP[i].track++] = lines[i].p2;
    }

    // Check which line every point lines on
    for (int i = 0; i < points_track; i++){
        for (int j = 0; j < line_track; j++){
            double lineX = lines[j].p2.x - lines[j].p1.x;
            double lineY = lines[j].p2.y - lines[j].p1.y;

            double pointX = lines[j].p2.x - points[i].x;
            double pointY = lines[j].p2.y - points[i].y;

            // Check if the point is on the line by using cross product
            bool pointOnLine = fabs(crossProduct(pointX, pointY, lineX, lineY)) < 0.5f;

            // Check if the point lines within the line
            bool pointInLine = points[i].x >= fmin(lines[j].p1.x, lines[j].p2.x) &&
                               points[i].x <= fmax(lines[j].p1.x, lines[j].p2.x) &&
                               points[i].y >= fmin(lines[j].p1.y, lines[j].p2.y) &&
                               points[i].y <= fmax(lines[j].p1.y, lines[j].p2.y);


            // Check if the point is an endpoint, if yes then dont add it
            bool isEndpoint = equalPoints(lineP[j].points[0], points[i]) || equalPoints(lineP[j].points[1], points[i]);

            // printf("Point inserted: (%f, %f)\n", points[i].x, points[i].y);
            // printf("Booleans: %d, %d, %d\n", pointOnLine, pointInLine, !isEndpoint);
            // The angles are almost equal and the point lies on the SEGMENTED LINE
            if (pointOnLine && pointInLine && !isEndpoint){
                lineP[j].points[lineP[j].track++] = points[i];

                //printf("Track: %d Size: %d\n", lineP[j].track, lineP[j].size);

                if (lineP[j].track >= lineP[j].size){
                    lineP[j].size *= 2;
                    resize((void**)&lineP[j].points, lineP[j].size, sizeof(SDL_FPoint));
                }
            }
        }
    }

    // Split every line
    for (int i = 0; i < line_track; i++){
        // if (lineP[i].track <= 2){
        //     addPoint(hashtable, lineP[i].points[0], lineP[i].points[1], pointInfo, adjInfo);
        //     addPoint(hashtable, lineP[i].points[1], lineP[i].points[0], pointInfo, adjInfo);
        // }

        for (int j = 0; j < lineP[i].track; j++){
            //printf("Point: (%f, %f)\n", lineP[i].points[j].x, lineP[i].points[j].y);
        }

        splitLine(hashtable, lineP[i].points, lineP[i].track, lines[i], newLines, new_track, 
                  new_size, pointInfo, adjInfo);
    }

    // Free the allocated memory
    for (int i = 0; i < line_track; i++){
        free(lineP[i].points);
    }

    free(lineP);
}

// Find all the free points starting from curP (assume curP is already a free endpoint)
void freeDFS(KeyNode* hashtable[], SDL_FPoint curP, bool* invalidLines, Line* lines, int line_track, 
             TypeInfo pointInfo){
    
    ValueNode* hashSearch = hashtableSearch(hashtable, &curP, pointInfo);

    // No neighbours (idk how this would happen but maybe it happens)
    if (!hashSearch){
        return;
    }
    AdjLines adjList = *(AdjLines*)hashSearch->value;

    int numValid = 0;
    int validIndexLine, validIndexAdj;
    // See if only one line is valid, then its a free point
    for (int i = 0; i < adjList.track; i++){
        int lineIndex = getLineIndex(lines, line_track, (Line){curP, adjList.adjPoints[i]});

        if (!invalidLines[lineIndex]){
            numValid++;
            validIndexLine = lineIndex;
            validIndexAdj = i;
        }
    }

    if (numValid == 1){
        invalidLines[validIndexLine] = true;
        freeDFS(hashtable, adjList.adjPoints[validIndexAdj], invalidLines, lines, line_track, pointInfo);
    }

}

// Finds all the connected nodes of an adjacency list starting from some point curP, and sets them to invalid
void connectedDFS(KeyNode* hashtable[], SDL_FPoint curP, bool* invalidLines, Line* lines, int line_track, TypeInfo pointInfo){
            
    // Extract the adjacent points from the hashtable
    ValueNode* hashSearch = hashtableSearch(hashtable, &curP, pointInfo);
    if (!hashSearch){
        return;
    }
    AdjLines adjList = *(AdjLines*)hashSearch->value;

    // Base case - all adjacent nodes have been processed
    for (int i = 0; i < adjList.track; i++){

        int sameLineIndex = getLineIndex(lines, line_track, (Line){curP, adjList.adjPoints[i]});

        // Line is still valid, there is still something left to check
        if (!invalidLines[sameLineIndex]){
            // Set the line to invalid
            invalidLines[sameLineIndex] = true;
            
            // Recursively check all nodes
            connectedDFS(hashtable, adjList.adjPoints[i], invalidLines, lines, line_track, pointInfo);
        }
    }
}

// Return bool so we can check if there is a polygon that actually enclosed our point or not
bool findPolygon(KeyNode* hashtable[], Line* lines, int line_track, int line_size, SDL_FPoint clicked, Line* polygon, 
                 int* polygon_track, int* polygon_size, TypeInfo pointInfo, TypeInfo adjInfo){
    
    // Find the region that encloses our point

    // Stores all the lines that we "remove"
    bool* invalidLines = malloc(line_size * sizeof(bool));
    for (int i = 0; i < line_track; i++){
        invalidLines[i] = false;
    }

    // Remove all lines that have atleast one free endpoint
    for (int i = 0; i < line_track; i++){
        freeDFS(hashtable, lines[i].p1, invalidLines, lines, line_track, pointInfo);
        freeDFS(hashtable, lines[i].p2, invalidLines, lines, line_track, pointInfo);
    }

    while (true){
        // Get the closest line segment
        float minDist = INFINITY;
        int minIndex = -1;
        for (int i = 0; i < line_track; i++){
            // If its a line we "removed", then skip it
            if (invalidLines[i]){
                continue;
            }

            SDL_FPoint p1 = lines[i].p1;
            SDL_FPoint p2 = lines[i].p2;

            // Get the line vector
            double dx = p1.x - p2.x;
            double dy = p1.y - p2.y;
            double length_squared = dx*dx + dy*dy;
            double distance;

            if (length_squared == 0.0) {
                // A and B are the same point
                distance = SDL_sqrt((clicked.x - p1.x) * (clicked.x - p1.x) + (clicked.y - p1.y) * (clicked.y - p1.y));
            }
            else {
                // Compute the projection scalar t
                double t = ((clicked.x - p2.x)*dx + (clicked.y - p2.y)*dy) / length_squared;

                // Clamp t to the range [0,1]
                if (t < 0.0) t = 0.0;
                else if (t > 1.0) t = 1.0;

                // Compute the projection point on the segment
                SDL_FPoint projection = {p2.x + t * dx, p2.y + t * dy};

                distance = SDL_sqrt((clicked.x - projection.x) * (clicked.x - projection.x) + 
                                    (clicked.y - projection.y) * (clicked.y - projection.y));
            }

            if (distance < minDist){
                minDist = distance;
                minIndex = i;
            }
        }

        if (debug){
            printf("Closest line segment: Point 1: (%f, %f), Point 2: (%f, %f)\n", lines[minIndex].p1.x, lines[minIndex].p1.y, 
                                                                                   lines[minIndex].p2.x, lines[minIndex].p2.y);
        }

        if (minIndex == -1){
            free(polygon);
            free(invalidLines);
            if (debug){
                printf("Couldn't find a region enclosing this point!\n");
            }
            return false;
        }

        // Traverse neighbouring lines until we find the initial line

        // Connect point - the point that is connected with the previous line
        // Next point - the point that will be connected to the next line
        SDL_FPoint connectPoint, nextPoint, startPoint;
        connectPoint = lines[minIndex].p1;
        nextPoint = lines[minIndex].p2;

        // Check if the point lies to the right of the vector, if yes then swap the points
        // (we need the point to lie to the left of the VECTOR)

        float cross = (nextPoint.x - connectPoint.x) * (clicked.y - connectPoint.y) - (nextPoint.y - connectPoint.y) * 
                      (clicked.x - connectPoint.x);

        // Clicked is on the right side of the vector (positive)
        if (cross > 0){
            // Swap endpoints (reverse the vector so that our point now lies to the left of it)
            SDL_FPoint swap = connectPoint;
            connectPoint = nextPoint;
            nextPoint = swap;
        }

        startPoint = connectPoint;

        // Unclear what should be filled
        if (cross == 0){
            return false;
        }

        *polygon_track = 0;

        polygon[*polygon_track].p1 = connectPoint;
        polygon[(*polygon_track)++].p2 = nextPoint;

        // Loop until we form a closed loop
        while (!(equalPoints(startPoint, nextPoint))){
            if (*polygon_track >= *polygon_size){
                (*polygon_size) *= 2;
                resize((void**)&polygon, *polygon_size, sizeof(Line));
            }

            // Prevent infinite looping (idk why it happens lowk)
            if (*polygon_size > 1600){
                // this print statement is very silly hehe
                printf("Polygon size limit hit! (1600 sides)");
                free(polygon);
                free(invalidLines);
                return false;
            }

            // Calculate the angle between next point and current point
            double dist_x = nextPoint.x - connectPoint.x;
            double dist_y = nextPoint.y - connectPoint.y;
            double lineAngle = atan2(dist_y, dist_x) * 180.0 / PI;
            lineAngle = fmod(lineAngle + 360.0, 360.0);

            // Loop through all adjacent line segments
            ValueNode* hashSearch = hashtableSearch(hashtable, &nextPoint, pointInfo);
            AdjLines adjList = *(AdjLines*)hashSearch->value;

            angleIndex* angleArray = malloc(adjList.track * sizeof(angleIndex));
            int angle_track = 0;

            for (int i = 0; i < adjList.track; i++){
                // Same point as the current lines other endpoint (we dont wanna check the same line!!!)
                if (equalPoints(adjList.adjPoints[i], connectPoint)){
                    continue;
                }

                // Check if its in invalid line and if yes, dont add it to angleArray (super duper important!!!)
                int lineIndex = getLineIndex(lines, line_track, (Line){adjList.adjPoints[i], nextPoint});

                if (invalidLines[lineIndex]){
                    continue;
                }

                // Calculate each angle
                double dist_x = adjList.adjPoints[i].x - nextPoint.x;
                double dist_y = adjList.adjPoints[i].y - nextPoint.y;
                double angle = atan2(dist_y, dist_x) * 180.0 / PI;
                angle = fmod(angle + 360.0, 360.0);

                // Straight line disconnected by some points
                double diff;
                if (fabs(angle - lineAngle) < 1e-2){
                    diff = 0;
                }
                else {
                    diff = 360 - fmod((angle - lineAngle + 360.0), 360.0);       // Gets the counter clockwise distance
                }

                angleArray[angle_track].angle = diff;
                angleArray[angle_track++].index = i;
            }

            qsort(angleArray, angle_track, sizeof(angleIndex), qsortCompareAngle);

            int angleIndex = angle_track-1;
            int bestIndex = 0;

            // Case where angle == 180 -> straight line, choose that
            // Find the biggest angle that is less than 180 (if you turn more than 180 youre going clockwise now)
            while (angleIndex > 0 && angleArray[angleIndex].angle > 180.0){
                angleIndex--;
            }

            bestIndex = angleArray[angleIndex].index;

            free(angleArray);

            // Move onto a new line by defining new endpoints
            connectPoint = nextPoint;
            nextPoint = adjList.adjPoints[bestIndex];

            //printf("Point 1: (%f, %f), Point 2: (%f, %f)\n", connectPoint.x, connectPoint.y, nextPoint.x, nextPoint.y);

            polygon[*polygon_track].p1 = connectPoint;
            polygon[(*polygon_track)++].p2 = nextPoint;
        }

        // Check if our point is inside the polygon
        bool inside = false;
        for (int i = 0; i < *polygon_track; i++){
            // Silly formula
            bool intersect = ((polygon[i].p2.y > clicked.y) != (polygon[i].p1.y > clicked.y)) &&
                            (clicked.x < (polygon[i].p1.x - polygon[i].p2.x) * (clicked.y - polygon[i].p2.y) / 
                                    (polygon[i].p1.y - polygon[i].p2.y) + polygon[i].p2.x);
            if (intersect)
                inside = !inside;
        }

        // Move on to the next part if the clicked part is inside the polygon
        if (inside){
            break;
        }

        // If not inside, remove all connected parts to the lines we just processed

        // Helper function to set all connected parts to invalid
        connectedDFS(hashtable, polygon[0].p1, invalidLines, lines, line_track, pointInfo);

        // Check if all lines are invalid (either already processed or pruned at the start) and break out of the function
        // if so
        bool allInvalid = true;
        for (int i = 0; i < line_track; i++){
            // Break immediately if we find a line that has not been processed
            if (!invalidLines[i]){
                allInvalid = false;
                break;
            }
        }

        if (allInvalid){
            free(polygon);
            free(invalidLines);
            if (debug){
                printf("Couldn't find a region enclosing this point!\n");
            }
            return false;
        }
    }

    return true;
}

// Fills the section that was clicked on - returns if a region was actually clicked or not
bool fillLines(SDL_Renderer* renderer, KeyNode* hashtable[], SDL_Color color, SDL_FPoint clicked, Line* lines, 
               int line_track, int line_size, TypeInfo pointInfo, TypeInfo adjInfo, float zoom, SDL_FPoint top_left){

    // Find the region that encloses our point

    Line* polygon = malloc(100 * sizeof(Line));
    int polygon_track = 0; int polygon_size = 100;

    // Find the closest enclosing polygon
    bool isEnclosed = findPolygon(hashtable, lines, line_track, line_size, clicked, polygon, &polygon_track, 
                                  &polygon_size, pointInfo, adjInfo);

    // Didnt find an enclosing polygon
    if (!isEnclosed){
        return false;
    }

    // Found the enclosing lines for a region, render using scanlines
    // Convert to screen cordinates
    for (int i = 0; i < polygon_track; i++){
        polygon[i].p1 = world_to_screen(polygon[i].p1, zoom, top_left);
        polygon[i].p2 = world_to_screen(polygon[i].p2, zoom, top_left);
    }

    float minY = 0;
    float maxY = WINDOW_HEIGHT;

    // Find the maximum and minimum X values of the polygon so we know where to draw the scanline
    float minX = INFINITY;
    float maxX = -INFINITY;

    for (int i = 0; i < polygon_track; i++){
        if (polygon[i].p1.x < minX) minX = polygon[i].p1.x;
        if (polygon[i].p2.x < minX) minX = polygon[i].p2.x;
        if (polygon[i].p1.x > maxX) maxX = polygon[i].p1.x;
        if (polygon[i].p2.x > maxX) maxX = polygon[i].p2.x; 
    }

    // Horizontal scanline starting point
    Line scanLine;
    scanLine.p1 = (SDL_FPoint){.x = minX, .y = minY};
    scanLine.p2 = (SDL_FPoint){.x = maxX, .y = minY};

    // Adjust the minY and maxY to fit inside the polygon
    minY++;
    maxY--;

    // Draw loop
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
    while (true){
        scanLine.p1.y += 1;
        scanLine.p2.y += 1;

        SDL_FPoint* intersects = malloc(polygon_track * sizeof(SDL_FPoint));
        int track = 0;
        int size = polygon_track;

        SDL_FPoint lineIntersects[1];

        // Find intersections with the scanline
        for (int i = 0; i < polygon_track; i++){
            bool hasIntersected = line_line_intersect(scanLine, polygon[i], lineIntersects);

            // Check if there was an intersection at all
            if (hasIntersected){
                intersects[track++] = lineIntersects[0];
            }
        }

        // Sort the intersections
        qsort(intersects, track, sizeof(SDL_FPoint), qsortCompPoint);

        // Adjust the intersection points so it doesn't draw over the lines
        for (int i = 0; i < track; i++){
            if (i % 2 == 0){
                intersects[i].x++;
            }
            else {
                intersects[i].x--;
            }
        }

        // Draw the lines by drawing between every pair of intersection points
        for (int i = 0; i < track; i += 2){
            Line intersectLine = (Line){.p1.x = intersects[i].x, .p1.y = intersects[i].y, .p2.x = intersects[(i+1) % track].x,
                                        .p2.y = intersects[(i+1) % track].y};

            drawClippedLine(renderer, intersectLine);
        }

        // End condition
        if (scanLine.p1.y >= maxY){
            free(intersects);
            break;
        }
    }

    free(polygon);

    return true;
}

// Unfill the region that was clicked on
void unfillLines(KeyNode* hashtable[], Line* lines, int line_track, int line_size, SDL_FPoint clicked, fillRegions* regions,
                 int* fill_track, TypeInfo pointInfo, TypeInfo adjInfo){

    // Find the region that we just clicked
    Line* clicked_p = malloc(100 * sizeof(Line));
    int clicked_track = 0; int clicked_size = 100;

    bool isEnclosedClick = findPolygon(hashtable, lines, line_track, line_size, clicked, clicked_p, &clicked_track, 
                                  &clicked_size, pointInfo, adjInfo);

    // No region was clicked
    if (!isEnclosedClick){
        return;
    }

    // Sort the array so they match
    qsort(clicked_p, clicked_track, sizeof(Line), qsortCompareLines);

    Line* region_p = malloc(100 * sizeof(Line));
    int region_track = 0; int region_size = 100;

    // Compare to all the regions 
    for (int i = 0; i < *fill_track; i++){
        bool isEnclosedRegion = findPolygon(hashtable, lines, line_track, line_size, regions[i].clicked, region_p,
                                            &region_track, &region_size, pointInfo, adjInfo);

        // The points in the fillRegion array are guaranteed to be enclosed by some region, but check anyway
        if (!isEnclosedRegion){
            continue;
        }

        // Sort this array too
        qsort(region_p, region_track, sizeof(Line), qsortCompareLines);

        // Compare
        // Same number of points?
        if (region_track != clicked_track){
            continue;
        }

        // Same cordinate points?
        bool sameCords = true;
        for (int j = 0; j < region_track; j++){
            // Not the same cordinates
            if (!(equalPoints(clicked_p[j].p1, region_p[j].p1) && equalPoints(clicked_p[j].p2, region_p[j].p2) &&
                !(equalPoints(clicked_p[j].p2, region_p[j].p1) && equalPoints(clicked_p[j].p1, region_p[j].p2)))){

                sameCords = false;
                break;
            }
        }

        // We have found the region, but there may be more points drawing this region so dont break yet
        if (sameCords){
            delete_arr(regions, i, *fill_track, sizeof(fillRegions));
            (*fill_track)--;

            // Because we deleted an element the indexes have shifted once to the left
            i--;
        }
    }
}

void drawThickLine(SDL_Renderer* renderer, Line line, int thickness){
    for (int j = 0; j <= thickness; j++){
        // Get a vector of the line and normalise it
        double distX = line.p1.x - line.p2.x;
        double distY = line.p1.y - line.p2.y;

        double dist = SDL_sqrt(distX * distX + distY * distY);

        distX /= dist;
        distY /= dist;

        // Rotate the vector
        double newX = -distY;
        double newY = distX;

        Line drawLine = line;
        // Adjust the line
        drawLine.p1.x += (j - thickness / 2) * newX;
        drawLine.p1.y += (j - thickness / 2) * newY;

        drawLine.p2.x += (j - thickness / 2) * newX;
        drawLine.p2.y += (j - thickness / 2) * newY;

        drawClippedLine(renderer, drawLine);
    }
}

int SDL_main(int argc, char* argv[]){

    SDL_Renderer* renderer;
    SDL_Window* window;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer("Construction simulation", WINDOW_WIDTH, WINDOW_HEIGHT, 0, &window, &renderer);
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    TTF_Init();

    // Holds properties
    char* texts[] = {"Currently using (D): ", "Radius: ", "Radius locked? (R): ", "Move around? (T): ", 
                     "Zoom amount: ", "Line type (H): ", "Fill colors? (K): ", "Color type (D): "};
    int size = sizeof(texts) / sizeof(texts[0]);

    // Holds the answers
    char** texts2 = (char**)malloc(size * sizeof(char*));

    SDL_Color textColor = {60, 127, 255, 255};
    SDL_Color textColor2 = {60, 255, 127, 255};
    SDL_Color fileColor = {255, 0, 0, 255};

    Tex* hold_texs = (Tex*)malloc(size * sizeof(Tex));

    Properties temp = 0;

    TTF_Font* font = TTF_OpenFont("C:\\Users\\user\\C_files\\tuff_dir\\Roboto_Condensed-Black.ttf", 20);
    TTF_Font* big_font =  TTF_OpenFont("C:\\Users\\user\\C_files\\tuff_dir\\Roboto_Condensed-Black.ttf", 40);

    for (int i = 0; i < size; i++){
        char* text = texts[i];

        hold_texs[i].pos = (SDL_FRect){.x = 0, .y = 0};

        hold_texs[i].tex = renderTexture(renderer, font, text, &hold_texs[i].pos, textColor);

        hold_texs[i].props = temp++;
    }

    SDL_Event event;
    bool run = true;

    // Store points
    SDL_FPoint* points = (SDL_FPoint*)malloc(100 * sizeof(SDL_FPoint));
    int points_track = 0;
    int points_size = 100;

    // For zooming and moving around
    SDL_FPoint* screen_points = (SDL_FPoint*)malloc(100 * sizeof(SDL_FPoint));

    // Stores areas that have been filled
    fillRegions* fillClicks = (fillRegions*)malloc(100 * sizeof(fillRegions));
    int fill_track = 0;
    int fill_size = 100;

    Line* lines = (Line*)malloc(100 * sizeof(Line));
    int line_track = 0;
    int line_size = 100;

    Line* screen_lines = (Line*)malloc(100 * sizeof(Line));

    Line* noIntersections = (Line*)malloc(100 * sizeof(Line));
    int noInt_track = 0; int noInt_size = 100;

    Circle* circles = (Circle*)malloc(100 * sizeof(Circle));
    int circle_track = 0;
    int circle_size = 100;

    Circle* screen_circles = (Circle*)malloc(100 * sizeof(Circle));

    State cur_state = DRAW_LINE;

    Dir rotate_dir = CLOCKWISE;

    bool shift_down = false;
    bool helper_line = false;
    bool fillColors = false;
    bool compassFinished = false;
    bool lockRadius = false;

    bool moveAround = false;
    bool firstClick = true;
    bool leftMouseDown = false;

    bool loadClicked = false;
    bool saveClicked = false;
    bool saveAsClicked = false;
    bool deleteClicked = false;

    int lineThickness = 1;
    int clicks_amount = 0;
    int help_line_alpha = 127;          // How light the helper lines appear
    int frame_count = 0;
    int fillIndex = 0;

    float radius = 10;                  // The locked radius amount
    float fps = 60.0f;

    float zoom = 1;
    double elapsed_time;

    char filename[256] = {0};

    SDL_FPoint prev_mouse = {0, 0};
    SDL_FPoint change = {0, 0};
    SDL_FPoint top_left = {0, 0};               // The top left corner of the viewport

    // Create a hashtable to store what points are adjacent to which ones
    KeyNode* hashtable[MAX_SIZE] = {NULL};

    TypeInfo pointInfo, adjInfo;
    initType(&pointInfo, sizeof(SDL_FPoint), cmp_point, hash_point, free_point, print_point);
    initType(&adjInfo, sizeof(AdjLines), cmp_adj, hash_adj, free_adj, print_adj);

    SDL_Color line_color = {.r = 0, .g = 255, .b = 0, .a = -1};

    // Initialize stuff
    SDL_Color colors[9] = {
        {255, 0, 0, 255},   // Red
        {180, 255, 180, 255},   // Green
        {0, 0, 255, 255},   // Blue
        {255, 255, 0, 255}, // Yellow
        {255, 165, 0, 255}, // Orange
        {128, 0, 128, 255}, // Purple
        {0, 255, 255, 255}, // Cyan
        {255, 192, 203, 255}, // Pink
        {255, 255, 255, 255}, // White
    };

    clock_t last_time = clock();
    clock_t start_time = clock();
    clock_t current_time;

    circles[circle_track].radius = radius;
    
    SDL_FRect saveAsRect;
    SDL_FRect saveRect;
    SDL_FRect loadRect;
    SDL_FRect deleteRect;

    while (run){

        SDL_SetRenderDrawColor(renderer, 160, 160, 160, 255);
        SDL_RenderClear(renderer);

        while (SDL_PollEvent(&event)){
            switch (event.type){
                case SDL_EVENT_QUIT:
                    run = false;
                    break;
                case SDL_EVENT_MOUSE_BUTTON_DOWN:
                    if (event.button.button == SDL_BUTTON_LEFT){
                        leftMouseDown = true;
                        firstClick = true;

                        if (moveAround){
                            break;
                        }

                        if (fillColors){
                            SDL_FPoint clicked = screen_to_world((SDL_FPoint){event.button.x, event.button.y}, zoom, top_left);
                            SDL_Color cur_color = colors[fillIndex];

                            // We need to remove the filled region
                            if (shift_down){
                                unfillLines(hashtable, noIntersections, noInt_track, noInt_size, clicked, fillClicks, 
                                            &fill_track, pointInfo, adjInfo);
                            }
                            // We need to add a filled region
                            else {
                                // Make sure the user actually clicked a region
                                Line* enclosingP = malloc(100 * sizeof(Line));
                                int enclosing_track = 0; int enclosing_size = 100;

                                bool isEnclosed = findPolygon(hashtable, noIntersections, noInt_track, noInt_size, clicked, 
                                                              enclosingP, &enclosing_track, &enclosing_size, pointInfo, adjInfo);

                                if (isEnclosed){
                                    fillClicks[fill_track++] = (fillRegions){.clicked = clicked, .color = cur_color};

                                    if (fill_track >= fill_size){
                                        fill_size *= 2;
                                        resize((void**)&fillClicks, fill_size, sizeof(fillRegions));
                                    }
                                }
                                else {
                                    printf("No region clicked!\n");
                                }

                                if (debug){
                                    printf("\n\n");
                                    fillLines(renderer, hashtable, fillClicks[fill_track-1].color, fillClicks[fill_track-1].clicked, 
                                              noIntersections, noInt_track, noInt_size, pointInfo, adjInfo, zoom, top_left);
                                }
                            }
                            break;
                        }

                        // Check if the "load" button was clicked
                        if (SDL_PointInRectFloat(&(SDL_FPoint){event.button.x, event.button.y}, &loadRect)){
                            loadClicked = true;
                        }

                        // Check if the "save" button was clicked
                        if (SDL_PointInRectFloat(&(SDL_FPoint){event.button.x, event.button.y}, &saveRect)){
                            saveClicked = true;
                        }

                        // Check if the "save as" button was clicked
                        if (SDL_PointInRectFloat(&(SDL_FPoint){event.button.x, event.button.y}, &saveAsRect)){
                            saveAsClicked = true;
                        }

                        // Check if the "delete" button was clicked
                        if (SDL_PointInRectFloat(&(SDL_FPoint){event.button.x, event.button.y}, &deleteRect)){
                            deleteClicked = true;
                        }

                        if (saveAsClicked){
                            break;
                        }

						// Increment the clicks amount so that I can use it when closing the tab that comes up
                        // when you press load, save or delete
                        if (loadClicked || saveClicked || deleteClicked){
                            clicks_amount %= 2;
                            clicks_amount++;

                            break;
                        }

                        switch (cur_state){
                            case NUM_ELEMENTS:
                                break;
                            case DRAW_LINE:
                                // Add support to remove lines
                                // Logic is identical to snapping a line to a point while drawing
                                if (shift_down){
                                    float minAngle = 361.0;
                                    int minIndex;

                                    for (int i = 0; i < line_track; i++){
                                        float point_dist_x = screen_lines[i].p1.x - screen_lines[i].p2.x;
                                        float point_dist_y = screen_lines[i].p1.y - screen_lines[i].p2.y;

                                        float point_dist = SDL_sqrt(point_dist_x * point_dist_x + point_dist_y * point_dist_y);

                                        float point_angle = atan2(point_dist_y, point_dist_x) * 180.0 / PI;
                                        point_angle = fmod(point_angle + 360.0, 360.0);

                                        float mouse_dist_x = screen_lines[i].p1.x - event.button.x;
                                        float mouse_dist_y = screen_lines[i].p1.y - event.button.y;

                                        float mouse_dist = SDL_sqrt(mouse_dist_x * mouse_dist_x + mouse_dist_y * mouse_dist_y);

                                        float mouse_angle = atan2(mouse_dist_y, mouse_dist_x) * 180.0 / PI;
                                        mouse_angle = fmod(mouse_angle + 360.0, 360.0);

                                        if (fabs(mouse_angle - point_angle) < minAngle && mouse_dist < point_dist){
                                            minAngle = fabs(mouse_angle - point_angle);
                                            minIndex = i;
                                        }
                                    }

                                    if (minAngle < 10){
                                        delete_arr((void*)lines, minIndex, line_track, sizeof(Line));
                                        line_track--;
                                    }

                                    break;
                                }
                                clicks_amount %= 2;
                                clicks_amount++;

                                // Just check for intersections (and increment), because you dont need to set the end point.
                                // It is already done in the loop that runs every frame.
                                // Also add this node to the hashtable that stores adjacency values. Either append to an
                                // already existing value, or create a new one.
                                if (clicks_amount == 2){
                                    points_track++;

                                    SDL_FPoint intersect;

                                    // Check for line line intersections
                                    for (int i = 0; i < line_track; i++){

                                        bool isIntersected = line_line_intersect(lines[line_track], lines[i], &intersect);

                                        if (isIntersected){
                                            points[points_track++] = intersect;

                                            if (points_track >= points_size){
                                                points_size *= 2;
                                                resize((void**)(&points), points_size, sizeof(SDL_FPoint));
                                                resize((void**)(&screen_points), points_size, sizeof(SDL_FPoint));
                                            }
                                        }
                                    }

                                    // Check for circle - line intersections
                                    SDL_FPoint intersections[2];
                                    int track = 0;
                                    for (int i = 0; i < circle_track; i++){
                                        circle_line_intersect(circles[i], lines[line_track], intersections, &track);

                                        for (int j = 0; j < track; j++){

                                            points[points_track++] = intersections[j];

                                            if (points_track >= points_size){
                                                points_size *= 2;
                                                resize((void**)(&points), points_size, sizeof(SDL_FPoint));
                                                resize((void**)(&screen_points), points_size, sizeof(SDL_FPoint));
                                            }
                                        }
                                    }

                                    line_track++;

                                    if (line_track >= line_size){
                                        line_size *= 2;

                                        resize((void**)&lines, line_size, sizeof(Line));
                                        resize((void**)&screen_lines, line_size, sizeof(Line));
                                    }

                                    if (noInt_track >= noInt_size){
                                        noInt_size *= 2;

                                        resize((void**)&noIntersections, noInt_size, sizeof(Line));
                                    }
                                }

                                // Set the start point, but make sure to account for zoom and moving around.
                                if (clicks_amount == 1){
                                    SDL_FPoint cur_clicked = screen_to_world((SDL_FPoint){event.button.x, event.button.y}, 
                                                                             zoom, top_left);

                                    lines[line_track].p1 = nearest_point(points, points_track, cur_clicked.x, cur_clicked.y,
                                                                         zoom);

                                    points[points_track++] = lines[line_track].p1;

                                    if (points_track >= points_size){
                                        points_size *= 2;
                                        resize((void**)&points, points_size, sizeof(SDL_FPoint));
                                        resize((void**)&screen_points, points_size, sizeof(SDL_FPoint));
                                    }
                                }
                                break;
                            case USE_COMPASS:
                                // Label must not be followed by a declaration!!! so ill just add this
                                line_track = line_track;

                                // The clicked point in absolute cordinates
                                SDL_FPoint cur_clicked = screen_to_world((SDL_FPoint){event.button.x, event.button.y}, zoom,
                                                                         top_left);

                                // Apply the nearest point function always
                                SDL_FPoint real_point = nearest_point(points, points_track, cur_clicked.x, cur_clicked.y
                                                                      ,zoom);

                                // Add support for removing circle segments 
                                // Same as the line removal, but different algorithm to account for start angle,
                                // end angle and radius

                                // There is one small problem with this: When clicking on the overlapping part
                                // of two circle segments, the algorithm removes the one which was placed last.
                                // This is just a consequence of how the code is written tho, and you shouldn't have
                                // to place two overlapping circle segments during a normal construction anyway.
                                double minRadius = INFINITY;
                                int minIndex;
                                if (shift_down){
                                    for (int i = 0; i < circle_track; i++){

                                        // Distance between the absolute mid point and absolute mouse position
                                        double dist_x = circles[i].mid_p.x - cur_clicked.x;
                                        double dist_y = circles[i].mid_p.y - cur_clicked.y;

                                        double mouse_radius = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);
                                        double angle = 180.0 + atan2(dist_y, dist_x) * 180.0 / PI;

                                        double dist_radius = fabs(mouse_radius - circles[i].radius);
                                        bool angleBetween;
                                        // Full circle
                                        if (fabs(circles[i].start_angle - circles[i].end_angle) < 1e-2){
                                            angleBetween = true;
                                        }
                                        else {
                                            angleBetween = isBetween(circles[i].start_angle, circles[i].end_angle, 
                                                                     angle, circles[i].rotate);
                                        }

                                        if (angleBetween && dist_radius <= minRadius){
                                            minIndex = i;
                                            minRadius = dist_radius;
                                        }
                                    }

                                    if (minRadius < 30){
                                        delete_arr(circles, minIndex, circle_track, sizeof(Circle));
                                        delete_arr(screen_circles, minIndex, circle_track, sizeof(Circle));
                                        circle_track--;
                                    }

                                    break;
                                }

                                clicks_amount++;
                                // Place mid point
                                if (clicks_amount == 1){
                                    circles[circle_track].mid_p = real_point;

                                    screen_circles[circle_track] = circles[circle_track];
                                    break;
                                }
                                
                                float dist_x = circles[circle_track].mid_p.x - real_point.x;
                                float dist_y = circles[circle_track].mid_p.y - real_point.y;

                                // Set the beginning of the arc
                                if (clicks_amount == 3){
                                    compassFinished = false;
                                    double angle = 180.0 + atan2(dist_y, dist_x) * 180.0 / PI;

                                    circles[circle_track].start_angle = angle;
                                    screen_circles[circle_track].start_angle = angle;
                                }

                                // Copy the previous circles things
                                if (clicks_amount > 3){
                                    circle_track++;
                                    if (circle_track >= circle_size){
                                        circle_size *= 2;
                                        resize((void**)(&circles), circle_size, sizeof(Circle));
                                        resize((void**)(&screen_circles), circle_size, sizeof(Circle));
                                    }

                                    compassFinished = false;

                                    circles[circle_track] = circles[circle_track-1];
                                    screen_circles[circle_track] = screen_circles[circle_track-1];

                                    double angle = 180.0 + atan2(dist_y, dist_x) * 180.0 / PI;

                                    circles[circle_track].start_angle = angle;
                                    screen_circles[circle_track].start_angle = angle;
                                }
                        }
                    }
                    else if (event.button.button == SDL_BUTTON_RIGHT){
                        SDL_FPoint cur_clicked = screen_to_world((SDL_FPoint){event.button.x, event.button.y}, zoom, 
                                                                 top_left);

                        // Delete a point if shift is down
                        if (shift_down){
                            SDL_FPoint nearest = nearest_point(points, points_track, cur_clicked.x, cur_clicked.y, zoom);

                            for (int i = 0; i < points_track; i++){
                                if (equalPoints(points[i], nearest)){
                                    delete_arr((void*)points, i, points_track, sizeof(SDL_FPoint));
                                    points_track--;

                                    break;
                                }
                            }

                            break;
                        }
                        points[points_track++] = cur_clicked;

                        if (points_track >= points_size){
                            points_size *= 2;
                            resize((void**)(&points), points_size, sizeof(SDL_FPoint));
                            resize((void**)(&screen_points), points_size, sizeof(SDL_FPoint));
                        }
                    }
                    break;
                case SDL_EVENT_MOUSE_BUTTON_UP:
                    // End the arc
                    if (event.button.button == SDL_BUTTON_LEFT){
                        leftMouseDown = false;
                        if (clicks_amount < 3){
                            break;
                        }

                        compassFinished = true;

                        SDL_FPoint intersections[2];
                        int track = 0;

                        // Check intersections with any lines
                        for (int i = 0; i < line_track; i++){
                            circle_line_intersect(circles[circle_track], lines[i], intersections, &track);

                            for (int j = 0; j < track; j++){
                                points[points_track++] = intersections[j];

                                if (points_track >= points_size){
                                    points_size *= 2;
                                    resize((void**)(&points), points_size, sizeof(SDL_FPoint));
                                    resize((void**)(&screen_points), points_size, sizeof(SDL_FPoint));
                                }
                            }
                        }

                        // Check intersections with other circle segments
                        for (int i = 0; i < circle_track; i++){
                            circle_segment_intersection(circles[circle_track], circles[i], intersections, &track);

                            for (int j = 0; j < track; j++){
                                points[points_track++] = intersections[j];

                                if (points_track >= points_size){
                                    points_size *= 2;
                                    resize((void**)(&points), points_size, sizeof(SDL_FPoint));
                                    resize((void**)(&screen_points), points_size, sizeof(SDL_FPoint));
                                }
                            }
                        }
                    }
                    break;
                case SDL_EVENT_KEY_DOWN:
                    switch (event.key.key){
                        case SDLK_D:
                            if (fillColors){
                                fillIndex++;
                                fillIndex %= sizeof(colors) / sizeof(colors[0]);
                            }
                            else {
                                // wow!!! so cool! it changes states in two lines! it may be confusing to u tho :)
                                // Uses the way enums are represented (as integers) to do some cool stuff

                                cur_state++;
                                cur_state %= NUM_ELEMENTS;

                                clicks_amount = 0;
                            }
                            break;
                        case SDLK_LSHIFT:
                            // Used to remove points and lines
                            shift_down = true;
                            break;
                        case SDLK_H:
                            // Used to switch to a more lighter line
                            helper_line = !helper_line;
                            break;
                        case SDLK_RETURN:
                            // Used to end a compass draw event
                            clicks_amount = 0;
                            circle_track++;
                            break;
                        case SDLK_P:
                            // Used to change the direction of the line that is drawn around the circle
                            rotate_dir = rotate_dir == CLOCKWISE ? ANTI_CLOCKWISE : CLOCKWISE;
                            break;
                        case SDLK_R:
                            // Used to set the radius at a specific value
                            lockRadius = !lockRadius;
                            radius = circles[circle_track].radius;
                            break;
                        case SDLK_T:
                            moveAround = !moveAround;
                            break;
                        case SDLK_K:
                            fillColors = !fillColors;
                            break;
                        case SDLK_J:
                            printf("The hashtable: \n\n");
                            hashtablePrint(hashtable);
                            break;
                    }
                    break;
                case SDL_EVENT_KEY_UP:
                    if (event.key.key == SDLK_LSHIFT){
                        shift_down = false;
                    }
                    break;
                case SDL_EVENT_MOUSE_WHEEL:
                    // Get mouse state
                    SDL_PumpEvents();

                    float x, y;
                    SDL_GetMouseState(&x, &y);

                    SDL_FPoint mouseBeforeZoom = (SDL_FPoint){x, y};
                    mouseBeforeZoom = screen_to_world(mouseBeforeZoom, zoom, top_left);

                    zoom -= event.wheel.y * zoom / 10.0;

                    if (zoom < 1e-5){
                        zoom = 1e-5;
                    }

                    SDL_FPoint mouseAfterZoom = (SDL_FPoint){x, y};
                    mouseAfterZoom = screen_to_world(mouseAfterZoom, zoom, top_left);

                    top_left.x += mouseBeforeZoom.x - mouseAfterZoom.x;
                    top_left.y += mouseBeforeZoom.y - mouseAfterZoom.y;
                    break;
            }
        }

        // Get mouse state
        SDL_PumpEvents();

        float x, y;
        SDL_GetMouseState(&x, &y);

        // Remove nearly identical points
        for (int i = 0; i < points_track; i++){
            for (int j = i + 1; j < points_track; j++){
                while (j < points_track && fabs(points[i].x - points[j].x) < 1e-2 && fabs(points[i].y - points[j].y) < 1e-2){
                    delete_arr((void*)points, j, points_track, sizeof(SDL_FPoint));
                    points_track--;
                }
            }
        }

        // Split the lines so there is no intersections (required for fill tool)
        splitAllLines(hashtable, points, points_track, &noIntersections, &noInt_track, &noInt_size, 
                      lines, screen_lines, line_track, line_size, pointInfo, adjInfo);

        // Delete superrrr small lines
        for (int i = 0; i < noInt_track; i++){
            while (i < noInt_track && fabs(noIntersections[i].p1.x - noIntersections[i].p2.x) < 1e-2 && 
                                      fabs(noIntersections[i].p1.y - noIntersections[i].p2.y) < 1e-2){

                delete_arr(noIntersections, i, noInt_track, sizeof(Line));
                removePoint(hashtable, noIntersections[i].p1, noIntersections[i].p2, pointInfo, adjInfo);
                removePoint(hashtable, noIntersections[i].p1, noIntersections[i].p2, pointInfo, adjInfo);

                noInt_track--;
            }
        }

        // Render the filled regions
        for (int i = 0; i < fill_track; i++){
            if (debug){
                break;
            }

            fillLines(renderer, hashtable, fillClicks[i].color, fillClicks[i].clicked, noIntersections, noInt_track, 
                      noInt_size, pointInfo, adjInfo, zoom, top_left);

            SDL_SetRenderDrawColor(renderer, 0, 127, 255, 255);
            drawPoint(renderer, world_to_screen(fillClicks[i].clicked, zoom, top_left), 5);
        }

        // Draw the points
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        for (int i = 0; i < points_track; i++){
            screen_points[i] = world_to_screen((SDL_FPoint){points[i].x, points[i].y}, zoom, top_left);

            // Only draw if its inside the boundary
            if ((screen_points[i].x + 5 < 0 || screen_points[i].x - 5 > WINDOW_WIDTH) || 
                (screen_points[i].y + 5 < 0 || screen_points[i].y - 5 > WINDOW_HEIGHT)){
                    
                continue;
            }

            drawPoint(renderer, screen_points[i], 5);
        }

        // Draw the lines
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
        for (int i = 0; i < line_track; i++){
            screen_lines[i].p1 = world_to_screen(lines[i].p1, zoom, top_left);
            screen_lines[i].p2 = world_to_screen(lines[i].p2, zoom, top_left);

            int linesAmount;
            
            if (lines[i].helper_line){
                linesAmount = lineThickness;
            }
            else {
                linesAmount = lineThickness * 2;
            }

            drawThickLine(renderer, screen_lines[i], linesAmount);
        }

        // Draw the segmented circles
        for (int i = 0; i < circle_track; i++){
            screen_circles[i].mid_p = world_to_screen(circles[i].mid_p, zoom, top_left);
            screen_circles[i].radius = circles[i].radius / zoom;

            SDL_SetRenderDrawColor(renderer, 0, 127, 255, circles[i].helper_line ? help_line_alpha : 255);

            Circle new_circle = screen_circles[i];
            new_circle.radius -= 2;

            for (int i = 0; i < 4; i++){
                drawClippedSegmentCircles(renderer, new_circle);

                new_circle.radius += 1;
            }
        }

        // Render all the UI and stuff
        char str[100];

        // FPS string 
        strcpy(str, "FPS: ");
        SDL_FRect dstRect = {.x = 10, .y = 10};

        SDL_Texture* tex = renderTexture(renderer, font, str, &dstRect, textColor);
        SDL_RenderTexture(renderer, tex, NULL, &dstRect);

        SDL_DestroyTexture(tex);

        // The FPS amount
        sprintf(str, "%.2f", fps);
        SDL_FRect new_dstRect = {.x = dstRect.x + dstRect.w + 10, .y = 10};

        tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);
        SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

        SDL_DestroyTexture(tex);

        // Folder icon
        char* bmp_path;
        SDL_asprintf(&bmp_path, "%s%s", SDL_GetBasePath(), "saveIcon.bmp");
        SDL_Surface* surface = SDL_LoadBMP(bmp_path);

        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, surface);

        int width = surface->w / 8.0;
        int height = surface->h / 8.0;

        SDL_DestroySurface(surface);

        // Place it on the screen
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 10, .y = 10,
                             .w = width, .h = height};

        SDL_RenderTexture(renderer, texture, NULL, &dstRect);

        // Save as icon
        saveAsRect = dstRect;

        // Save icon - same texture, different place
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 10, .y = height + 20, .w = width, .h = height};

        saveRect = dstRect;

        SDL_RenderTexture(renderer, texture, NULL, &dstRect);

        // Load icon

        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 10, .y = 2 * height + 30, .w = width, .h = height};

        loadRect = dstRect;

        SDL_RenderTexture(renderer, texture, NULL, &dstRect);

        // Delete icon

        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 10, .y = 3 * height + 40, .w = width, .h = height};

        deleteRect = dstRect;

        SDL_RenderTexture(renderer, texture, NULL, &dstRect);

        SDL_DestroyTexture(texture);

        // Load text
        strcpy(str, "Load ->");
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 140, .y = 2 * height + 30};

        tex = renderTexture(renderer, big_font, str, &dstRect, textColor2);
        SDL_RenderTexture(renderer, tex, NULL, &dstRect);

        SDL_DestroyTexture(tex);

        // Save text
        strcpy(str, "Save ->");
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 140, .y = height + 20};

        tex = renderTexture(renderer, big_font, str, &dstRect, textColor2);
        SDL_RenderTexture(renderer, tex, NULL, &dstRect);

        SDL_DestroyTexture(tex);

        // Save as text
        strcpy(str, "Save as ->");
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - width - 190, .y = 10};

        tex = renderTexture(renderer, big_font, str, &dstRect, textColor2);
        SDL_RenderTexture(renderer, tex, NULL, &dstRect);

        SDL_DestroyTexture(tex);

        // Delete text
        strcpy(str, "Delete ->");
        dstRect = (SDL_FRect){.x = WINDOW_WIDTH - 215, .y = 3 * height + 40};

        tex = renderTexture(renderer, big_font, str, &dstRect, textColor2);
        SDL_RenderTexture(renderer, tex, NULL, &dstRect);

        SDL_DestroyTexture(tex);

        SDL_free(bmp_path);

        // Display frames
        current_time = clock();

        elapsed_time = (double)(current_time - start_time) / CLOCKS_PER_SEC;

        if (elapsed_time >= 1.0f) {
            fps = (float)(frame_count) / (elapsed_time);

            frame_count = 0;
            start_time = clock();  // Reset the timer
        }

        frame_count++;

        // Render properties
        for (int i = 0; i < size; i++){
            SDL_FRect dstRect;

            if (i == 0){
                hold_texs[i].pos.x = 10;
                hold_texs[i].pos.y = 50;

                dstRect = (SDL_FRect){10, 50, hold_texs[i].pos.w, hold_texs[i].pos.h};
                SDL_RenderTexture(renderer, hold_texs[i].tex, NULL, &dstRect);
            }
            else {
                int pos_y = hold_texs[i-1].pos.y + hold_texs[i-1].pos.h + 30;
                hold_texs[i].pos.y = pos_y;
                hold_texs[i].pos.x = 10;

                dstRect = (SDL_FRect){10, pos_y, hold_texs[i].pos.w, hold_texs[i].pos.h};
                SDL_RenderTexture(renderer, hold_texs[i].tex, NULL, &dstRect);
            }

            // Render the answers
            char str[100];
            SDL_Texture* tex;

            SDL_FRect new_dstRect = hold_texs[i].pos;
            new_dstRect.x += new_dstRect.w + 10;

            switch (hold_texs[i].props){
                case USING:
                    if (cur_state == DRAW_LINE){
                        strcpy(str, "Straight edge");

                        tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                        SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                        SDL_DestroyTexture(tex);
                    }
                    if (cur_state == USE_COMPASS){
                        strcpy(str, "Compass");

                        tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                        SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                        SDL_DestroyTexture(tex);
                    }
                    break;
                case RADIUS:
                    sprintf(str, "%.2f", circles[circle_track].radius);

                    tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case LOCKED:
                    if (lockRadius){
                        strcpy(str, "true");
                    }
                    else {
                        strcpy(str, "false");
                    }

                    tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case MOVING:
                    if (moveAround){
                        strcpy(str, "true");
                    }
                    else {
                        strcpy(str, "false");
                    }

                    tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case ZOOM:
                    sprintf(str, "%f", zoom);

                    tex = renderTexture(renderer, font , str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case HELP_LINE:
                    if (helper_line){
                        strcpy(str, "Helper line");
                    }
                    else {
                        strcpy(str, "Solid line");
                    }

                    tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case FILL:
                    if (fillColors){
                        strcpy(str, "true");
                    }
                    else {
                        strcpy(str, "false");
                    }

                    tex = renderTexture(renderer, font, str, &new_dstRect, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &new_dstRect);

                    SDL_DestroyTexture(tex);
                    break;
                case CLR_TYPE:
                    // Render a rectangle based on fillIndex
                    if (fillColors){
                        new_dstRect.h = 30;
                        new_dstRect.w = 30;

                        SDL_SetRenderDrawColor(renderer, colors[fillIndex].r, colors[fillIndex].g, colors[fillIndex].b, 255);

                        SDL_RenderFillRect(renderer, &new_dstRect);
                    }
            }
        }

        // Save the game data to a variable at first, then choose where to write it
        if (saveClicked || saveAsClicked){
            SaveData fileData;

            // Allocate memory for the arrays
            fileData.points.arr = (SDL_FPoint*)malloc(points_track * sizeof(SDL_FPoint));
            fileData.screen_points.arr = (SDL_FPoint*)malloc(points_track * sizeof(SDL_FPoint));

            fileData.lines.arr = (Line*)malloc(line_track * sizeof(Line));
            fileData.screen_lines.arr = (Line*)malloc(line_track * sizeof(Line));

            fileData.circles.arr = (Circle*)malloc(circle_track * sizeof(Circle));
            fileData.screen_circles.arr = (Circle*)malloc(circle_track * sizeof(Circle));

            fileData.fillRegions.arr = (fillRegions*)malloc(fill_track * sizeof(fillRegions));

            // Copy the data
            // Points
            memcpy(fileData.points.arr, points, points_track * sizeof(SDL_FPoint));
            fileData.points.track = points_track;
            fileData.points.size = points_size;

            memcpy(fileData.screen_points.arr, screen_points, points_track * sizeof(SDL_FPoint));
            fileData.screen_points.track = points_track;
            fileData.screen_points.size = points_size;

            // Lines
            memcpy(fileData.lines.arr, lines, line_track * sizeof(Line));
            fileData.lines.track = line_track;
            fileData.lines.size = line_size;

            memcpy(fileData.screen_lines.arr, screen_lines, line_track * sizeof(Line));
            fileData.screen_lines.track = line_track;
            fileData.screen_lines.size = line_size;

            // Circles
            memcpy(fileData.circles.arr, circles, circle_track * sizeof(Circle));
            fileData.circles.track = circle_track;
            fileData.circles.size = circle_size;

            memcpy(fileData.screen_circles.arr, screen_circles, circle_track * sizeof(Circle));
            fileData.screen_circles.track = circle_track;
            fileData.screen_circles.size = circle_size;

            // Fill regions
            memcpy(fileData.fillRegions.arr, fillClicks, fill_track * sizeof(fillRegions));
            fileData.fillRegions.track = fill_track;
            fileData.fillRegions.size = fill_size;

            // Screen data
            fileData.zoom = zoom;
            fileData.top_left = top_left;

            SDL_Event newEvent;

            if (saveAsClicked){
                // Get a file name
                SDL_StartTextInput(window);
                int index = 0;
                while (SDL_PollEvent(&newEvent)) {
                    // Get the characters the user is typing
                    if (newEvent.type == SDL_EVENT_TEXT_INPUT){
                        if (strlen(filename) + strlen(newEvent.text.text) < sizeof(filename))
                            strcat(filename, newEvent.text.text);
                    }
                    // Handle writing to the file when user is done naming it
                    else if (newEvent.type == SDL_EVENT_KEY_DOWN && newEvent.key.key == SDLK_RETURN){
                        SDL_StopTextInput(window);

                        // Copy the name into the entry
                        strcpy(fileData.name, filename);

                        // Append it to the end of the file
                        saveAFile(&fileData);

                        saveAsClicked = false;
                        break;
                    }
                    // Remove characters
                    else if (newEvent.type == SDL_EVENT_KEY_DOWN && newEvent.key.key == SDLK_BACKSPACE){
                        if (strlen(filename) > 0){
                            filename[strlen(filename)-1] = '\0';
                        }
                    }
                }

                // Render the name youre writing
                if (strlen(filename) > 0){
                    SDL_FRect pos = {.x = 10, .y = WINDOW_HEIGHT - 100};
                    SDL_Texture* tex = renderTexture(renderer, font, filename, &pos, textColor2);

                    SDL_RenderTexture(renderer, tex, NULL, &pos);

                    SDL_DestroyTexture(tex);
                }
            }
            // Write to the file you are CURRENTLY ON!!!
            else {
                saveClicked = false;
                
                strcpy(fileData.name, filename);

                // Save to the right file
                saveToFile(&fileData, filename);

                clicks_amount = 0;
                leftMouseDown = false;
            }
        }

        // Display all the file data
        if (loadClicked){
            SaveData* dataArray[MAX_FILES];
            int count = loadAllSaveFiles(dataArray, MAX_FILES);

            if (count == 0){
                loadClicked = false;
                clicks_amount = 0;
                continue;
            }

            // Render the data array
            
            // Store all the positions of each entry (on the screen)
            SDL_FRect* positions = (SDL_FRect*)malloc(count * sizeof(SDL_FRect));
            positions[0] = (SDL_FRect){.x = 0, .y = 4 * height + 40, .w = 0, .h = 0};     // Adjust initial Y to fit correctly
            for (int i = 0; i < count; i++){
                char* string = dataArray[i]->name;

                SDL_FRect pos;
                SDL_Texture* tex = renderTexture(renderer, font, string, &pos, fileColor);

                pos.x = WINDOW_WIDTH - pos.w - 10;
                if (i == 0){
                    pos.y = positions[i].y + 10;
                }
                else {
                    pos.y = positions[i-1].y + positions[i-1].h + 10;
                }

                positions[i] = pos;

                SDL_RenderTexture(renderer, tex, NULL, &pos);

                SDL_DestroyTexture(tex);
            }

            if (leftMouseDown){
                // Check if the load button was clicked again
                if (SDL_PointInRectFloat(&(SDL_FPoint){x,y}, &loadRect) && clicks_amount >= 2){

                    loadClicked = false;
                    clicks_amount = 0;
                }
                for (int i = 0; i < count; i++){
                    SDL_FRect pos = positions[i];

                    // Check if the file entry was clicked
                    if (SDL_PointInRectFloat(&(SDL_FPoint){x,y}, &pos)){
                        loadClicked = false;
                        // Load all the data into the game arrays

                        // Load all the data into the game arrays from dataArray[i]
                        SaveData* selected = dataArray[i];

                        // Save the NAME!!!!
                        strcpy(filename, selected->name);

                        // Points
                        memcpy(points, selected->points.arr, selected->points.track * sizeof(SDL_FPoint));
                        points_track = selected->points.track;
                        points_size = selected->points.size;

                        memcpy(screen_points, selected->screen_points.arr, selected->screen_points.track * sizeof(SDL_FPoint));

                        // Lines
                        memcpy(lines, selected->lines.arr, selected->lines.track * sizeof(Line));
                        line_track = selected->lines.track;
                        line_size = selected->lines.size;

                        memcpy(screen_lines, selected->screen_lines.arr, selected->screen_lines.track * sizeof(Line));

                        // Circles
                        memcpy(circles, selected->circles.arr, selected->circles.track * sizeof(Circle));
                        circle_track = selected->circles.track;
                        circle_size = selected->circles.size;

                        memcpy(screen_circles, selected->screen_circles.arr, selected->screen_circles.track * sizeof(Circle));

                        // Fill regions
                        memcpy(fillClicks, selected->fillRegions.arr, selected->fillRegions.track * sizeof(fillRegions));
                        fill_track = selected->fillRegions.track;
                        fill_size = selected->fillRegions.size;

                        // Screen data
                        zoom = selected->zoom;
                        top_left = selected->top_left;
                    }
                }
            }
            
            free(positions);
            freeData(dataArray, count);
        }

        if (deleteClicked){
            SaveData* dataArray[MAX_FILES];
            int count = loadAllSaveFiles(dataArray, MAX_FILES);

            if (count == 0){
                deleteClicked = false;
                continue;
            }

            // Render the data array
            
            // Store all the positions of each entry (on the screen)
            SDL_FRect* positions = (SDL_FRect*)malloc(count * sizeof(SDL_FRect));
            positions[0] = (SDL_FRect){.x = 0, .y = 4 * height + 40, .w = 0, .h = 0};     // Adjust initial Y to fit correctly
            for (int i = 0; i < count; i++){
                char* string = dataArray[i]->name;

                SDL_FRect pos;
                SDL_Texture* tex = renderTexture(renderer, font, string, &pos, fileColor);

                pos.x = WINDOW_WIDTH - pos.w - 10;
                if (i == 0){
                    pos.y = positions[i].y + 10;
                }
                else {
                    pos.y = positions[i-1].y + positions[i-1].h + 10;
                }

                positions[i] = pos;

                SDL_RenderTexture(renderer, tex, NULL, &pos);

                SDL_DestroyTexture(tex);
            }

            if (leftMouseDown){
                // Check if the load button was clicked again
                if (SDL_PointInRectFloat(&(SDL_FPoint){x,y}, &deleteRect) && clicks_amount >= 2){

                    deleteClicked = false;
                    clicks_amount = 0;
                }
                for (int i = 0; i < count; i++){
                    SDL_FRect pos = positions[i];

                    // Check if the file entry was clicked
                    if (SDL_PointInRectFloat(&(SDL_FPoint){x,y}, &pos)){
                        deleteClicked = false;

                        deleteFile(dataArray[i]->name);

                        leftMouseDown = false;
                        break;
                    }
                }
            }
            
            free(positions);
            freeData(dataArray, count);
        }

        // After all the rendering is done, dont go further
        if (loadClicked || saveAsClicked || saveClicked || deleteClicked){
            SDL_RenderPresent(renderer);
            continue;
        }

        // The absolute cordinates in the world
        SDL_FPoint cur_clicked = screen_to_world((SDL_FPoint){x, y}, zoom, top_left);

        // Display the line or circle segment you want to delete
        if (shift_down){
            // Display the region we want to delete
            if (fillColors){
                SDL_Color white = {0, 0, 0, 100};
                fillLines(renderer, hashtable, white, cur_clicked, noIntersections, noInt_track, noInt_size, pointInfo, 
                          adjInfo, zoom, top_left);
            }
            // Display the point you want to delete aswell
            SDL_FPoint nearest = nearest_point(screen_points, points_track, x, y, 1);
            
            // Different color, since the point is already red
            SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
            if (nearest.x != x && nearest.y != y){
                drawPoint(renderer, nearest, 5);
            }

            // Render the adjacent points (as blue)
            SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
            SDL_FPoint searchP = screen_to_world(nearest, zoom, top_left);
            ValueNode* hashSearch = hashtableSearch(hashtable, &searchP, pointInfo);

            if (hashSearch){
                AdjLines adjList = *(AdjLines*)hashSearch->value;

                for (int i = 0; i < adjList.track; i++){
                    drawPoint(renderer, world_to_screen(adjList.adjPoints[i], zoom, top_left), 5);
                }
            }

            SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
            switch (cur_state){
                case DRAW_LINE: {
                    float minAngle = 361.0;
                    int minIndex;

                    for (int i = 0; i < line_track; i++){
                        float point_dist_x = screen_lines[i].p1.x - screen_lines[i].p2.x;
                        float point_dist_y = screen_lines[i].p1.y - screen_lines[i].p2.y;

                        float point_dist = SDL_sqrt(point_dist_x * point_dist_x + point_dist_y * point_dist_y);

                        float point_angle = atan2(point_dist_y, point_dist_x) * 180.0 / PI;
                        point_angle = fmod(point_angle + 360.0, 360.0);

                        float mouse_dist_x = screen_lines[i].p1.x - x;
                        float mouse_dist_y = screen_lines[i].p1.y - y;

                        float mouse_dist = SDL_sqrt(mouse_dist_x * mouse_dist_x + mouse_dist_y * mouse_dist_y);

                        float mouse_angle = atan2(mouse_dist_y, mouse_dist_x) * 180.0 / PI;
                        mouse_angle = fmod(mouse_angle + 360.0, 360.0);

                        if (fabs(mouse_angle - point_angle) < minAngle && mouse_dist < point_dist){
                            minAngle = fabs(mouse_angle - point_angle);
                            minIndex = i;
                        }
                    }

                    if (minAngle < 10){
                        Line cur = screen_lines[minIndex];

                        int linesAmount = lines[minIndex].helper_line ? lineThickness : lineThickness * 2;

                        drawThickLine(renderer, cur, linesAmount);
                    }
                    break;
                }
                case USE_COMPASS: {
                    double minRadius = INFINITY;
                    int minIndex;
                    if (shift_down){
                        for (int i = 0; i < circle_track; i++){

                            // Distance between the absolute mid point and absolute mouse position
                            double dist_x = circles[i].mid_p.x - cur_clicked.x;
                            double dist_y = circles[i].mid_p.y - cur_clicked.y;

                            double mouse_radius = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);
                            double angle = 180.0 + atan2(dist_y, dist_x) * 180.0 / PI;

                            double dist_radius = fabs(mouse_radius - circles[i].radius);
                            bool angleBetween;
                            // Full circle
                            if (fabs(circles[i].start_angle - circles[i].end_angle) < 1e-2){
                                angleBetween = true;
                            }
                            else {
                                angleBetween = isBetween(circles[i].start_angle, circles[i].end_angle, 
                                                            angle, circles[i].rotate);
                            }

                            if (angleBetween && dist_radius <= minRadius){
                                minIndex = i;
                                minRadius = dist_radius;
                            }
                        }

                        if (minRadius < 30){
                            Circle new_circle = screen_circles[minIndex];
                            new_circle.radius -= 2;

                            for (int i = 0; i < 4; i++){
                                drawClippedSegmentCircles(renderer, new_circle);

                                new_circle.radius += 1;
                            }
                        }

                        break;
                    }
                }
            }
        }

        if (moveAround && leftMouseDown){
            // The current mouse position
            SDL_FPoint cur_mouse = {x * zoom, y * zoom};

            // First click is special (everything will be relative to this)
            if (firstClick){
                firstClick = false;

                prev_mouse = cur_mouse;
            }

            // This is more like a vector but whatever vectors are a kind of point after all

            change.x = cur_mouse.x - prev_mouse.x;
            change.y = cur_mouse.y - prev_mouse.y;

            // Set current mouse pos to previous (in the next iteration it will be the previous)
            prev_mouse = cur_mouse;

            top_left.x -= change.x;
            top_left.y -= change.y;
        }

        // Draws the line as you change it
        if (cur_state == DRAW_LINE && clicks_amount == 1){
            SDL_FPoint clicked = lines[line_track].p1;
            SDL_FPoint draw, screen_draw;

            lines[line_track].helper_line = helper_line;

            screen_draw = nearest_point(screen_points, points_track, x, y, 1);      // No zoom

            draw = screen_to_world(screen_draw, zoom, top_left);

            SDL_SetRenderDrawColor(renderer, line_color.r, line_color.g, line_color.b, 255);

            // Snapping the line to the closest point
            float minAngle = INFINITY;
            SDL_FPoint minPoint;
            // Check for intersections with the changing line
            for (int i = 0; i < points_track; i++){
                if (points[i].x == clicked.x && points[i].y == clicked.y){
                    continue;
                }

                float point_dist_x = screen_lines[line_track].p1.x - screen_points[i].x;
                float point_dist_y = screen_lines[line_track].p1.y - screen_points[i].y;

                float point_dist = SDL_sqrt(point_dist_x * point_dist_x + point_dist_y * point_dist_y);

                float point_angle = atan2(point_dist_y, point_dist_x) * 180.0 / PI;

                float mouse_dist_x = screen_lines[line_track].p1.x - screen_draw.x;
                float mouse_dist_y = screen_lines[line_track].p1.y - screen_draw.y;

                float mouse_dist = SDL_sqrt(mouse_dist_x * mouse_dist_x + mouse_dist_y * mouse_dist_y);

                float mouse_angle = atan2(mouse_dist_y, mouse_dist_x) * 180.0 / PI;

                // Find the smallest angle, because there may be multiple angles that are less than 10
                // Also make sure that it doesn't lock onto the point unless the mouse is already past that point
                if (fabs(mouse_angle - point_angle) < minAngle && mouse_dist >= point_dist){
                    minPoint = points[i];
                    minAngle = fabs(mouse_angle - point_angle);
                }
            }

            if (minAngle < 10){
                // Calculate the vector and normalize it
                // Point from the first clicked on end to the lines another end
                float dist_x = minPoint.x - clicked.x;
                float dist_y = minPoint.y - clicked.y;

                float dist = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);

                dist_x /= dist;
                dist_y /= dist;

                // Scale the vector by the distance
                float new_dist_x = draw.x - clicked.x;
                float new_dist_y = draw.y - clicked.y;

                dist = SDL_sqrt(new_dist_x * new_dist_x + new_dist_y * new_dist_y);

                dist_x *= dist;
                dist_y *= dist;

                lines[line_track].p2 = (SDL_FPoint){clicked.x + dist_x, clicked.y + dist_y};
                points[points_track] = (SDL_FPoint){clicked.x + dist_x, clicked.y + dist_y};
            }
            else {
                lines[line_track].p2 = (SDL_FPoint){draw.x, draw.y};
                points[points_track] = (SDL_FPoint){draw.x, draw.y};
            }

            screen_lines[line_track].p1 = world_to_screen(lines[line_track].p1, zoom, top_left);
            screen_lines[line_track].p2 = world_to_screen(lines[line_track].p2, zoom, top_left);

            screen_points[points_track] = world_to_screen(points[points_track], zoom, top_left);

            int linesAmount = lines[line_track].helper_line ? lineThickness : lineThickness * 2;
            
            drawThickLine(renderer, screen_lines[line_track], linesAmount);
            // SDL_RenderLine(renderer, screen_lines[line_track].p1.x, screen_lines[line_track].p1.y,
            //                screen_lines[line_track].p2.x, screen_lines[line_track].p2.y);

            SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
            drawPoint(renderer, screen_points[points_track], 5);
        }

        // Draw the radius and circle as you change it
        if (cur_state == USE_COMPASS && clicks_amount > 0){
            circles[circle_track].helper_line = helper_line;
            circles[circle_track].rotate = rotate_dir;

            screen_circles[circle_track].helper_line = helper_line;
            screen_circles[circle_track].rotate = rotate_dir;

            SDL_SetRenderDrawColor(renderer, 0, 127, 255, helper_line ? help_line_alpha : 255);

            // The absolute position of the mouse in the world, applied with the nearest point function
            SDL_FPoint real_point = nearest_point(points, points_track, cur_clicked.x, cur_clicked.y, zoom);

            // Distance between the middle of the absolute circle and the mouse absolute position
            double dist_x = circles[circle_track].mid_p.x - real_point.x;
            double dist_y = circles[circle_track].mid_p.y - real_point.y;

            double dist = SDL_sqrt(dist_x * dist_x + dist_y * dist_y);

            SDL_FPoint screen_mid_p = world_to_screen(circles[circle_track].mid_p, zoom, top_left);

            // Re-calculate the mid point
            screen_circles[circle_track].mid_p = screen_mid_p;

            // Distance between the screen circles mid and the screen mouse point
            SDL_FPoint screen_nearest = world_to_screen(real_point, zoom, top_left);

            double screen_dist_x = screen_mid_p.x - screen_nearest.x;
            double screen_dist_y = screen_mid_p.y - screen_nearest.y;

            if (clicks_amount == 1){
                if (lockRadius){
                    circles[circle_track].radius = radius;
                    screen_circles[circle_track].radius = radius / zoom;        // Scaled down version of radius

                    // Normalize and scale
                    dist_x /= dist;
                    dist_y /= dist;

                    dist_x *= radius;
                    dist_y *= radius;

                    SDL_FPoint point2 = world_to_screen((SDL_FPoint){circles[circle_track].mid_p.x - dist_x, 
                                                        circles[circle_track].mid_p.y - dist_y}, zoom, top_left);

                    // Render in the direction of the mouse but only render certain radius
                    SDL_RenderLine(renderer, screen_mid_p.x, screen_mid_p.y, point2.x, point2.y);

                    drawCircle(renderer, screen_mid_p, screen_circles[circle_track].radius);
                }
                else {
                    circles[circle_track].radius = dist;
                    screen_circles[circle_track].radius = dist / zoom;      // Scaled down version of radius

                    SDL_RenderLine(renderer, screen_mid_p.x, screen_mid_p.y, screen_nearest.x, screen_nearest.y);

                    drawCircle(renderer, screen_mid_p, screen_circles[circle_track].radius);

                }
            }

            screen_circles[circle_track].radius = circles[circle_track].radius / zoom;

            if (clicks_amount >= 2){
                drawCircle(renderer, screen_mid_p, screen_circles[circle_track].radius);
            }

            if (clicks_amount >= 3){
                if (!compassFinished){
                    double end_angle = 180.0 + atan2(screen_dist_y, screen_dist_x) * 180.0 / PI;

                    circles[circle_track].end_angle = end_angle;
                    screen_circles[circle_track].end_angle = end_angle;
                }

                Circle new_circle = screen_circles[circle_track];
                new_circle.radius -= 2;

                for (int i = 0; i < 4; i++){
                    drawClippedSegmentCircles(renderer, new_circle);

                    new_circle.radius += 1;
                }
            }
        }

        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_QuitSubSystem(SDL_INIT_VIDEO);
    SDL_Quit();
    TTF_CloseFont(font);
    TTF_Quit();

    free(points);
    free(screen_points);
    free(lines);
    free(screen_lines);
    free(circles);
    free(screen_circles);
    free(hold_texs);
    free(texts2);

    return 0;

}
