#define VTX_CNT 88
#define EDGE_CNT 1117
int graf_vtx_cnt = VTX_CNT;
int graf_edge_cnt = EDGE_CNT;
int graf_edge_data[EDGE_CNT][3] = {
{ 1, 10, 1 },
{ 1, 16, 2 },
{ 1, 22, 3 },
{ 1, 5, 1 },
{ 1, 53, 1 },
{ 1, 60, 1 },
{ 1, 82, 1 },
{ 2, 16, 2 },
{ 2, 27, 1 },
{ 2, 39, 3 },
{ 2, 60, 1 },
{ 2, 7, 1 },
{ 2, 82, 1 },
{ 3, 16, 2 },
{ 3, 20, 2 },
{ 3, 22, 1 },
{ 3, 33, 2 },
{ 3, 37, 2 },
{ 3, 42, 2 },
{ 3, 55, 2 },
{ 3, 60, 3 },
{ 3, 64, 1 },
{ 3, 67, 2 },
{ 3, 7, 3 },
{ 3, 76, 2 },
{ 3, 79, 2 },
{ 3, 84, 1 },
{ 4, 20, 1 },
{ 4, 22, 2 },
{ 4, 39, 2 },
{ 5, 10, 1 },
{ 5, 1, 1 },
{ 5, 12, 1 },
{ 5, 16, 2 },
{ 5, 18, 1 },
{ 5, 21, 1 },
{ 5, 22, 3 },
{ 5, 24, 1 },
{ 5, 27, 1 },
{ 5, 28, 1 },
{ 5, 33, 2 },
{ 5, 34, 1 },
{ 5, 36, 1 },
{ 5, 37, 2 },
{ 5, 39, 3 },
{ 5, 42, 4 },
{ 5, 46, 2 },
{ 5, 48, 2 },
{ 5, 49, 2 },
{ 5, 52, 2 },
{ 5, 53, 1 },
{ 5, 54, 1 },
{ 5, 60, 1 },
{ 5, 6, 2 },
{ 5, 62, 1 },
{ 5, 63, 1 },
{ 5, 64, 3 },
{ 5, 65, 1 },
{ 5, 69, 1 },
{ 5, 7, 1 },
{ 5, 75, 2 },
{ 5, 76, 4 },
{ 5, 77, 4 },
{ 5, 79, 4 },
{ 5, 81, 2 },
{ 5, 82, 1 },
{ 5, 86, 1 },
{ 5, 87, 1 },
{ 5, 88, 2 },
{ 6, 16, 1 },
{ 6, 20, 1 },
{ 6, 22, 2 },
{ 6, 25, 2 },
{ 6, 27, 2 },
{ 6, 31, 1 },
{ 6, 33, 1 },
{ 6, 39, 2 },
{ 6, 42, 3 },
{ 6, 44, 1 },
{ 6, 50, 1 },
{ 6, 5, 2 },
{ 6, 53, 2 },
{ 6, 55, 3 },
{ 6, 60, 2 },
{ 6, 64, 3 },
{ 6, 67, 3 },
{ 6, 7, 2 },
{ 6, 76, 3 },
{ 6, 78, 3 },
{ 6, 82, 2 },
{ 6, 86, 2 },
{ 7, 10, 1 },
{ 7, 12, 1 },
{ 7, 14, 1 },
{ 7, 15, 1 },
{ 7, 17, 1 },
{ 7, 18, 1 },
{ 7, 2, 1 },
{ 7, 22, 3 },
{ 7, 27, 1 },
{ 7, 28, 1 },
{ 7, 32, 1 },
{ 7, 3, 3 },
{ 7, 33, 2 },
{ 7, 34, 2 },
{ 7, 35, 1 },
{ 7, 36, 1 },
{ 7, 37, 2 },
{ 7, 39, 3 },
{ 7, 40, 2 },
{ 7, 41, 1 },
{ 7, 42, 4 },
{ 7, 43, 1 },
{ 7, 46, 1 },
{ 7, 50, 1 },
{ 7, 52, 2 },
{ 7, 54, 1 },
{ 7, 56, 2 },
{ 7, 60, 1 },
{ 7, 6, 2 },
{ 7, 64, 3 },
{ 7, 66, 1 },
{ 7, 68, 1 },
{ 7, 69, 1 },
{ 7, 72, 1 },
{ 7, 74, 3 },
{ 7, 75, 2 },
{ 7, 76, 4 },
{ 7, 79, 4 },
{ 7, 8, 1 },
{ 7, 81, 2 },
{ 7, 82, 1 },
{ 7, 87, 1 },
{ 8, 16, 2 },
{ 8, 20, 2 },
{ 8, 33, 2 },
{ 8, 36, 1 },
{ 8, 39, 3 },
{ 8, 52, 2 },
{ 8, 60, 1 },
{ 8, 7, 1 },
{ 8, 81, 2 },
{ 8, 82, 1 },
{ 9, 22, 1 },
{ 9, 39, 1 },
{ 9, 55, 2 },
{ 9, 64, 1 },
{ 9, 80, 2 },
{ 10, 1, 1 },
{ 10, 16, 2 },
{ 10, 33, 2 },
{ 10, 37, 2 },
{ 10, 48, 2 },
{ 10, 50, 1 },
{ 10, 5, 1 },
{ 10, 60, 1 },
{ 10, 7, 1 },
{ 10, 81, 2 },
{ 11, 16, 1 },
{ 11, 20, 1 },
{ 11, 33, 1 },
{ 11, 60, 2 },
{ 12, 16, 2 },
{ 12, 39, 3 },
{ 12, 5, 1 },
{ 12, 60, 1 },
{ 12, 64, 3 },
{ 12, 7, 1 },
{ 12, 82, 1 },
{ 13, 22, 2 },
{ 13, 39, 2 },
{ 13, 55, 1 },
{ 13, 64, 2 },
{ 13, 73, 1 },
{ 13, 78, 1 },
{ 14, 16, 2 },
{ 14, 33, 2 },
{ 14, 50, 1 },
{ 14, 7, 1 },
{ 15, 16, 2 },
{ 15, 33, 2 },
{ 15, 60, 1 },
{ 15, 7, 1 },
{ 16, 10, 2 },
{ 16, 11, 1 },
{ 16, 1, 2 },
{ 16, 12, 2 },
{ 16, 13, 3 },
{ 16, 14, 2 },
{ 16, 15, 2 },
{ 16, 17, 2 },
{ 16, 18, 2 },
{ 16, 2, 2 },
{ 16, 22, 2 },
{ 16, 23, 1 },
{ 16, 24, 2 },
{ 16, 27, 2 },
{ 16, 28, 2 },
{ 16, 3, 2 },
{ 16, 32, 2 },
{ 16, 33, 2 },
{ 16, 35, 2 },
{ 16, 37, 1 },
{ 16, 39, 2 },
{ 16, 40, 1 },
{ 16, 42, 3 },
{ 16, 43, 2 },
{ 16, 45, 2 },
{ 16, 46, 2 },
{ 16, 49, 1 },
{ 16, 50, 2 },
{ 16, 51, 2 },
{ 16, 5, 2 },
{ 16, 52, 1 },
{ 16, 53, 2 },
{ 16, 55, 3 },
{ 16, 56, 1 },
{ 16, 57, 1 },
{ 16, 58, 2 },
{ 16, 59, 3 },
{ 16, 60, 2 },
{ 16, 6, 1 },
{ 16, 63, 2 },
{ 16, 64, 2 },
{ 16, 65, 2 },
{ 16, 67, 3 },
{ 16, 68, 2 },
{ 16, 69, 2 },
{ 16, 70, 2 },
{ 16, 7, 2 },
{ 16, 72, 2 },
{ 16, 73, 3 },
{ 16, 74, 2 },
{ 16, 76, 3 },
{ 16, 77, 3 },
{ 16, 78, 3 },
{ 16, 79, 3 },
{ 16, 80, 3 },
{ 16, 81, 1 },
{ 16, 8, 2 },
{ 16, 82, 2 },
{ 16, 83, 1 },
{ 16, 84, 2 },
{ 16, 85, 2 },
{ 16, 88, 1 },
{ 17, 16, 2 },
{ 17, 39, 3 },
{ 17, 50, 1 },
{ 17, 7, 1 },
{ 18, 16, 2 },
{ 18, 22, 3 },
{ 18, 27, 1 },
{ 18, 28, 1 },
{ 18, 39, 3 },
{ 18, 50, 1 },
{ 18, 5, 1 },
{ 18, 60, 1 },
{ 18, 64, 3 },
{ 18, 7, 1 },
{ 18, 81, 2 },
{ 18, 82, 1 },
{ 19, 33, 1 },
{ 20, 11, 1 },
{ 20, 25, 1 },
{ 20, 3, 1 },
{ 20, 33, 1 },
{ 20, 40, 1 },
{ 20, 4, 1 },
{ 20, 44, 1 },
{ 20, 47, 1 },
{ 20, 52, 1 },
{ 20, 56, 1 },
{ 20, 60, 1 },
{ 20, 6, 1 },
{ 20, 75, 1 },
{ 20, 8, 1 },
{ 20, 81, 1 },
{ 20, 83, 1 },
{ 20, 88, 1 },
{ 21, 22, 3 },
{ 21, 5, 1 },
{ 21, 60, 1 },
{ 21, 7, 1 },
{ 21, 82, 1 },
{ 22, 10, 3 },
{ 22, 1, 3 },
{ 22, 16, 2 },
{ 22, 18, 3 },
{ 22, 21, 3 },
{ 22, 24, 3 },
{ 22, 27, 3 },
{ 22, 29, 3 },
{ 22, 3, 1 },
{ 22, 32, 3 },
{ 22, 33, 2 },
{ 22, 34, 2 },
{ 22, 36, 2 },
{ 22, 37, 2 },
{ 22, 39, 1 },
{ 22, 4, 2 },
{ 22, 45, 3 },
{ 22, 48, 2 },
{ 22, 49, 2 },
{ 22, 50, 3 },
{ 22, 51, 3 },
{ 22, 52, 2 },
{ 22, 5, 3 },
{ 22, 53, 3 },
{ 22, 55, 2 },
{ 22, 56, 2 },
{ 22, 60, 3 },
{ 22, 6, 2 },
{ 22, 63, 3 },
{ 22, 64, 1 },
{ 22, 65, 3 },
{ 22, 67, 2 },
{ 22, 68, 3 },
{ 22, 69, 3 },
{ 22, 70, 1 },
{ 22, 7, 3 },
{ 22, 73, 2 },
{ 22, 74, 1 },
{ 22, 75, 2 },
{ 22, 77, 2 },
{ 22, 78, 2 },
{ 22, 79, 2 },
{ 22, 80, 2 },
{ 22, 81, 2 },
{ 22, 82, 1 },
{ 22, 83, 2 },
{ 22, 84, 1 },
{ 22, 86, 3 },
{ 22, 9, 1 },
{ 23, 39, 2 },
{ 24, 16, 2 },
{ 24, 22, 3 },
{ 24, 39, 3 },
{ 24, 50, 1 },
{ 24, 5, 1 },
{ 24, 60, 1 },
{ 24, 64, 3 },
{ 24, 7, 1 },
{ 24, 81, 2 },
{ 25, 20, 1 },
{ 25, 33, 1 },
{ 25, 39, 2 },
{ 25, 42, 3 },
{ 25, 6, 1 },
{ 25, 75, 3 },
{ 25, 76, 3 },
{ 26, 39, 3 },
{ 26, 60, 1 },
{ 26, 7, 1 },
{ 26, 82, 1 },
{ 27, 12, 1 },
{ 27, 16, 2 },
{ 27, 18, 1 },
{ 27, 2, 1 },
{ 27, 22, 3 },
{ 27, 32, 1 },
{ 27, 33, 2 },
{ 27, 34, 1 },
{ 27, 36, 1 },
{ 27, 37, 2 },
{ 27, 39, 3 },
{ 27, 41, 1 },
{ 27, 48, 2 },
{ 27, 50, 1 },
{ 27, 5, 1 },
{ 27, 52, 2 },
{ 27, 6, 2 },
{ 27, 63, 1 },
{ 27, 64, 3 },
{ 27, 65, 1 },
{ 27, 68, 1 },
{ 27, 69, 1 },
{ 27, 7, 1 },
{ 27, 81, 2 },
{ 27, 82, 1 },
{ 28, 16, 2 },
{ 28, 18, 1 },
{ 28, 32, 1 },
{ 28, 34, 1 },
{ 28, 48, 2 },
{ 28, 5, 1 },
{ 28, 63, 1 },
{ 28, 65, 1 },
{ 28, 7, 1 },
{ 28, 81, 2 },
{ 28, 86, 1 },
{ 29, 22, 3 },
{ 29, 60, 1 },
{ 29, 7, 1 },
{ 29, 81, 2 },
{ 29, 82, 1 },
{ 30, 16, 2 },
{ 30, 33, 2 },
{ 30, 60, 1 },
{ 30, 7, 1 },
{ 31, 20, 1 },
{ 31, 33, 1 },
{ 31, 6, 1 },
{ 31, 75, 1 },
{ 32, 16, 2 },
{ 32, 22, 3 },
{ 32, 27, 1 },
{ 32, 28, 1 },
{ 32, 39, 3 },
{ 32, 60, 1 },
{ 32, 82, 1 },
{ 33, 10, 2 },
{ 33, 11, 1 },
{ 33, 14, 2 },
{ 33, 15, 2 },
{ 33, 16, 1 },
{ 33, 19, 1 },
{ 33, 20, 1 },
{ 33, 22, 2 },
{ 33, 25, 2 },
{ 33, 27, 2 },
{ 33, 30, 2 },
{ 33, 31, 1 },
{ 33, 3, 2 },
{ 33, 34, 2 },
{ 33, 35, 2 },
{ 33, 36, 2 },
{ 33, 37, 1 },
{ 33, 39, 2 },
{ 33, 40, 1 },
{ 33, 42, 3 },
{ 33, 47, 1 },
{ 33, 50, 2 },
{ 33, 51, 2 },
{ 33, 5, 2 },
{ 33, 52, 1 },
{ 33, 53, 2 },
{ 33, 55, 3 },
{ 33, 56, 1 },
{ 33, 59, 3 },
{ 33, 60, 2 },
{ 33, 61, 2 },
{ 33, 64, 2 },
{ 33, 65, 2 },
{ 33, 69, 21 },
{ 33, 7, 2 },
{ 33, 75, 1 },
{ 33, 76, 3 },
{ 33, 79, 3 },
{ 33, 81, 2 },
{ 33, 8, 2 },
{ 33, 82, 2 },
{ 33, 83, 1 },
{ 34, 22, 2 },
{ 34, 27, 2 },
{ 34, 28, 2 },
{ 34, 33, 1 },
{ 34, 37, 1 },
{ 34, 42, 3 },
{ 34, 5, 2 },
{ 34, 60, 2 },
{ 34, 64, 2 },
{ 34, 7, 2 },
{ 34, 82, 2 },
{ 35, 16, 2 },
{ 35, 33, 2 },
{ 35, 60, 1 },
{ 35, 7, 1 },
{ 36, 16, 2 },
{ 36, 22, 3 },
{ 36, 27, 1 },
{ 36, 33, 2 },
{ 36, 39, 3 },
{ 36, 50, 1 },
{ 36, 5, 1 },
{ 36, 54, 1 },
{ 36, 7, 1 },
{ 36, 8, 1 },
{ 37, 10, 2 },
{ 37, 16, 1 },
{ 37, 20, 1 },
{ 37, 22, 2 },
{ 37, 27, 2 },
{ 37, 3, 2 },
{ 37, 33, 1 },
{ 37, 34, 2 },
{ 37, 39, 2 },
{ 37, 42, 3 },
{ 37, 48, 1 },
{ 37, 49, 1 },
{ 37, 5, 2 },
{ 37, 52, 1 },
{ 37, 55, 3 },
{ 37, 56, 1 },
{ 37, 60, 2 },
{ 37, 64, 2 },
{ 37, 67, 3 },
{ 37, 7, 2 },
{ 37, 76, 3 },
{ 37, 79, 3 },
{ 37, 81, 1 },
{ 37, 82, 2 },
{ 38, 52, 1 },
{ 38, 60, 1 },
{ 38, 82, 1 },
{ 39, 12, 3 },
{ 39, 13, 2 },
{ 39, 16, 2 },
{ 39, 17, 3 },
{ 39, 18, 3 },
{ 39, 2, 3 },
{ 39, 23, 2 },
{ 39, 24, 2 },
{ 39, 25, 1 },
{ 39, 26, 3 },
{ 39, 27, 3 },
{ 39, 3, 1 },
{ 39, 32, 3 },
{ 39, 33, 2 },
{ 39, 34, 2 },
{ 39, 36, 3 },
{ 39, 37, 2 },
{ 39, 40, 2 },
{ 39, 4, 2 },
{ 39, 42, 2 },
{ 39, 43, 3 },
{ 39, 44, 2 },
{ 39, 45, 3 },
{ 39, 47, 2 },
{ 39, 48, 2 },
{ 39, 50, 2 },
{ 39, 52, 2 },
{ 39, 5, 3 },
{ 39, 55, 2 },
{ 39, 56, 2 },
{ 39, 57, 2 },
{ 39, 58, 2 },
{ 39, 59, 2 },
{ 39, 60, 3 },
{ 39, 6, 2 },
{ 39, 63, 3 },
{ 39, 64, 1 },
{ 39, 65, 3 },
{ 39, 67, 2 },
{ 39, 69, 3 },
{ 39, 70, 1 },
{ 39, 7, 3 },
{ 39, 73, 2 },
{ 39, 74, 2 },
{ 39, 75, 2 },
{ 39, 76, 2 },
{ 39, 77, 2 },
{ 39, 78, 2 },
{ 39, 79, 2 },
{ 39, 80, 2 },
{ 39, 81, 2 },
{ 39, 82, 3 },
{ 39, 8, 3 },
{ 39, 83, 2 },
{ 39, 84, 2 },
{ 39, 88, 2 },
{ 39, 9, 1 },
{ 40, 16, 1 },
{ 40, 20, 1 },
{ 40, 33, 1 },
{ 40, 39, 2 },
{ 40, 64, 2 },
{ 40, 7, 2 },
{ 40, 81, 1 },
{ 41, 27, 1 },
{ 41, 60, 1 },
{ 41, 7, 1 },
{ 41, 82, 1 },
{ 41, 87, 1 },
{ 42, 16, 3 },
{ 42, 22, 2 },
{ 42, 25, 2 },
{ 42, 3, 2 },
{ 42, 33, 3 },
{ 42, 37, 3 },
{ 42, 48, 3 },
{ 42, 50, 4 },
{ 42, 52, 3 },
{ 42, 5, 4 },
{ 42, 55, 1 },
{ 42, 6, 3 },
{ 42, 64, 2 },
{ 42, 70, 2 },
{ 42, 73, 1 },
{ 42, 7, 4 },
{ 42, 74, 2 },
{ 42, 75, 3 },
{ 42, 77, 1 },
{ 42, 78, 1 },
{ 42, 81, 3 },
{ 42, 84, 2 },
{ 43, 16, 2 },
{ 43, 22, 3 },
{ 43, 39, 3 },
{ 43, 60, 1 },
{ 43, 64, 3 },
{ 43, 7, 1 },
{ 43, 81, 2 },
{ 43, 82, 1 },
{ 44, 20, 1 },
{ 44, 39, 2 },
{ 44, 6, 1 },
{ 45, 16, 2 },
{ 45, 22, 3 },
{ 45, 39, 3 },
{ 45, 60, 1 },
{ 45, 82, 1 },
{ 46, 16, 1 },
{ 46, 60, 2 },
{ 46, 7, 2 },
{ 47, 20, 1 },
{ 47, 33, 1 },
{ 47, 39, 2 },
{ 48, 10, 2 },
{ 48, 22, 2 },
{ 48, 27, 2 },
{ 48, 28, 2 },
{ 48, 37, 1 },
{ 48, 39, 2 },
{ 48, 42, 3 },
{ 48, 49, 1 },
{ 48, 5, 2 },
{ 48, 52, 1 },
{ 48, 53, 2 },
{ 48, 60, 2 },
{ 48, 64, 2 },
{ 48, 7, 2 },
{ 48, 77, 3 },
{ 48, 79, 3 },
{ 48, 81, 1 },
{ 48, 82, 1 },
{ 48, 85, 2 },
{ 49, 16, 1 },
{ 49, 22, 2 },
{ 49, 37, 1 },
{ 49, 48, 1 },
{ 49, 5, 2 },
{ 49, 60, 2 },
{ 49, 64, 2 },
{ 49, 81, 1 },
{ 50, 10, 1 },
{ 50, 14, 1 },
{ 50, 16, 2 },
{ 50, 17, 1 },
{ 50, 18, 1 },
{ 50, 22, 3 },
{ 50, 24, 1 },
{ 50, 27, 1 },
{ 50, 30, 1 },
{ 50, 33, 2 },
{ 50, 36, 1 },
{ 50, 37, 2 },
{ 50, 39, 3 },
{ 50, 42, 4 },
{ 50, 51, 1 },
{ 50, 52, 2 },
{ 50, 53, 1 },
{ 50, 54, 1 },
{ 50, 60, 1 },
{ 50, 61, 1 },
{ 50, 6, 2 },
{ 50, 63, 1 },
{ 50, 64, 3 },
{ 50, 65, 1 },
{ 50, 69, 1 },
{ 50, 7, 1 },
{ 50, 75, 2 },
{ 50, 76, 4 },
{ 50, 79, 4 },
{ 50, 81, 2 },
{ 50, 82, 1 },
{ 51, 16, 2 },
{ 51, 22, 3 },
{ 51, 33, 2 },
{ 51, 50, 1 },
{ 51, 52, 2 },
{ 51, 6, 2 },
{ 51, 64, 3 },
{ 51, 81, 2 },
{ 52, 16, 1 },
{ 52, 20, 1 },
{ 52, 22, 2 },
{ 52, 27, 2 },
{ 52, 33, 1 },
{ 52, 37, 1 },
{ 52, 38, 2 },
{ 52, 39, 2 },
{ 52, 42, 3 },
{ 52, 48, 1 },
{ 52, 50, 2 },
{ 52, 5, 2 },
{ 52, 53, 2 },
{ 52, 60, 2 },
{ 52, 64, 2 },
{ 52, 7, 2 },
{ 52, 81, 1 },
{ 52, 8, 2 },
{ 52, 82, 2 },
{ 53, 1, 1 },
{ 53, 16, 2 },
{ 53, 22, 3 },
{ 53, 33, 2 },
{ 53, 48, 2 },
{ 53, 50, 1 },
{ 53, 5, 1 },
{ 53, 81, 2 },
{ 54, 16, 2 },
{ 54, 36, 1 },
{ 54, 50, 1 },
{ 54, 5, 1 },
{ 54, 60, 1 },
{ 54, 82, 1 },
{ 55, 13, 1 },
{ 55, 16, 3 },
{ 55, 22, 2 },
{ 55, 3, 2 },
{ 55, 33, 3 },
{ 55, 37, 3 },
{ 55, 39, 2 },
{ 55, 42, 1 },
{ 55, 58, 1 },
{ 55, 59, 1 },
{ 55, 6, 3 },
{ 55, 64, 2 },
{ 55, 67, 1 },
{ 55, 74, 2 },
{ 55, 76, 1 },
{ 55, 79, 1 },
{ 55, 80, 1 },
{ 55, 9, 2 },
{ 56, 16, 1 },
{ 56, 22, 2 },
{ 56, 33, 1 },
{ 56, 37, 1 },
{ 56, 39, 2 },
{ 56, 64, 2 },
{ 56, 7, 2 },
{ 56, 81, 1 },
{ 57, 16, 1 },
{ 57, 22, 2 },
{ 57, 39, 2 },
{ 57, 64, 2 },
{ 57, 81, 1 },
{ 58, 16, 3 },
{ 58, 22, 2 },
{ 58, 39, 2 },
{ 58, 55, 1 },
{ 58, 73, 1 },
{ 58, 78, 1 },
{ 59, 16, 3 },
{ 59, 22, 2 },
{ 59, 33, 3 },
{ 59, 39, 2 },
{ 59, 55, 1 },
{ 59, 64, 2 },
{ 59, 73, 1 },
{ 59, 77, 1 },
{ 59, 78, 1 },
{ 60, 10, 1 },
{ 60, 1, 1 },
{ 60, 11, 2 },
{ 60, 12, 1 },
{ 60, 14, 1 },
{ 60, 16, 2 },
{ 60, 18, 1 },
{ 60, 2, 1 },
{ 60, 21, 1 },
{ 60, 22, 3 },
{ 60, 24, 1 },
{ 60, 29, 1 },
{ 60, 30, 1 },
{ 60, 32, 1 },
{ 60, 3, 3 },
{ 60, 33, 2 },
{ 60, 34, 1 },
{ 60, 35, 1 },
{ 60, 38, 1 },
{ 60, 39, 3 },
{ 60, 41, 1 },
{ 60, 43, 1 },
{ 60, 45, 1 },
{ 60, 46, 1 },
{ 60, 48, 2 },
{ 60, 49, 2 },
{ 60, 50, 1 },
{ 60, 5, 1 },
{ 60, 52, 2 },
{ 60, 54, 1 },
{ 60, 60, 1 },
{ 60, 6, 2 },
{ 60, 63, 1 },
{ 60, 64, 1 },
{ 60, 65, 1 },
{ 60, 68, 1 },
{ 60, 69, 1 },
{ 60, 7, 1 },
{ 60, 71, 1 },
{ 60, 72, 1 },
{ 60, 76, 4 },
{ 60, 8, 1 },
{ 60, 81, 2 },
{ 60, 86, 1 },
{ 61, 33, 2 },
{ 61, 50, 1 },
{ 61, 7, 1 },
{ 61, 81, 2 },
{ 62, 5, 1 },
{ 63, 16, 21 },
{ 63, 22, 3 },
{ 63, 27, 1 },
{ 63, 28, 1 },
{ 63, 39, 3 },
{ 63, 50, 1 },
{ 63, 5, 1 },
{ 63, 60, 1 },
{ 63, 64, 3 },
{ 63, 81, 2 },
{ 63, 82, 1 },
{ 63, 87, 1 },
{ 64, 12, 3 },
{ 64, 13, 2 },
{ 64, 16, 2 },
{ 64, 18, 3 },
{ 64, 22, 1 },
{ 64, 24, 3 },
{ 64, 25, 1 },
{ 64, 3, 1 },
{ 64, 33, 2 },
{ 64, 34, 3 },
{ 64, 37, 2 },
{ 64, 39, 1 },
{ 64, 40, 2 },
{ 64, 42, 2 },
{ 64, 48, 2 },
{ 64, 49, 2 },
{ 64, 50, 3 },
{ 64, 52, 2 },
{ 64, 5, 3 },
{ 64, 55, 2 },
{ 64, 56, 2 },
{ 64, 57, 2 },
{ 64, 58, 2 },
{ 64, 59, 2 },
{ 64, 60, 3 },
{ 64, 6, 2 },
{ 64, 63, 3 },
{ 64, 65, 3 },
{ 64, 67, 2 },
{ 64, 69, 3 },
{ 64, 70, 1 },
{ 64, 7, 3 },
{ 64, 73, 2 },
{ 64, 74, 1 },
{ 64, 75, 2 },
{ 64, 76, 2 },
{ 64, 77, 2 },
{ 64, 78, 2 },
{ 64, 79, 2 },
{ 64, 80, 2 },
{ 64, 81, 2 },
{ 64, 82, 3 },
{ 64, 83, 2 },
{ 64, 9, 1 },
{ 65, 16, 2 },
{ 65, 22, 3 },
{ 65, 27, 1 },
{ 65, 28, 1 },
{ 65, 33, 2 },
{ 65, 39, 3 },
{ 65, 50, 1 },
{ 65, 60, 1 },
{ 65, 64, 3 },
{ 65, 7, 1 },
{ 65, 82, 1 },
{ 66, 7, 1 },
{ 67, 16, 3 },
{ 67, 22, 2 },
{ 67, 3, 2 },
{ 67, 37, 3 },
{ 67, 39, 2 },
{ 67, 55, 1 },
{ 67, 6, 3 },
{ 67, 64, 2 },
{ 67, 70, 2 },
{ 67, 73, 1 },
{ 67, 74, 2 },
{ 67, 78, 1 },
{ 68, 16, 2 },
{ 68, 22, 3 },
{ 68, 27, 1 },
{ 68, 28, 1 },
{ 68, 60, 1 },
{ 68, 7, 1 },
{ 68, 82, 1 },
{ 69, 16, 2 },
{ 69, 22, 3 },
{ 69, 27, 1 },
{ 69, 33, 2 },
{ 69, 39, 3 },
{ 69, 50, 1 },
{ 69, 5, 1 },
{ 69, 60, 1 },
{ 69, 64, 3 },
{ 69, 7, 1 },
{ 69, 81, 2 },
{ 69, 82, 1 },
{ 70, 16, 3 },
{ 70, 22, 1 },
{ 70, 39, 1 },
{ 70, 42, 2 },
{ 70, 55, 2 },
{ 70, 64, 1 },
{ 70, 67, 2 },
{ 70, 76, 2 },
{ 70, 79, 2 },
{ 71, 5, 1 },
{ 71, 60, 1 },
{ 72, 16, 2 },
{ 72, 60, 1 },
{ 72, 7, 1 },
{ 72, 82, 1 },
{ 73, 13, 1 },
{ 73, 16, 3 },
{ 73, 22, 2 },
{ 73, 39, 2 },
{ 73, 42, 1 },
{ 73, 58, 1 },
{ 73, 59, 1 },
{ 73, 64, 2 },
{ 73, 67, 1 },
{ 73, 76, 1 },
{ 73, 79, 1 },
{ 74, 16, 2 },
{ 74, 22, 1 },
{ 74, 39, 1 },
{ 74, 42, 2 },
{ 74, 55, 2 },
{ 74, 64, 1 },
{ 74, 67, 2 },
{ 74, 7, 3 },
{ 74, 79, 2 },
{ 75, 16, 1 },
{ 75, 20, 1 },
{ 75, 22, 2 },
{ 75, 25, 2 },
{ 75, 31, 1 },
{ 75, 33, 1 },
{ 75, 39, 2 },
{ 75, 42, 3 },
{ 75, 50, 2 },
{ 75, 5, 2 },
{ 75, 60, 2 },
{ 75, 64, 2 },
{ 75, 7, 2 },
{ 75, 76, 2 },
{ 75, 81, 1 },
{ 75, 82, 2 },
{ 76, 16, 3 },
{ 76, 22, 2 },
{ 76, 3, 2 },
{ 76, 33, 3 },
{ 76, 37, 3 },
{ 76, 39, 2 },
{ 76, 50, 4 },
{ 76, 55, 1 },
{ 76, 6, 3 },
{ 76, 64, 2 },
{ 76, 70, 2 },
{ 76, 73, 1 },
{ 76, 7, 4 },
{ 76, 75, 3 },
{ 76, 77, 1 },
{ 76, 81, 3 },
{ 76, 84, 2 },
{ 77, 16, 3 },
{ 77, 22, 2 },
{ 77, 39, 2 },
{ 77, 42, 1 },
{ 77, 48, 3 },
{ 77, 5, 4 },
{ 77, 59, 1 },
{ 77, 64, 2 },
{ 77, 76, 1 },
{ 78, 13, 1 },
{ 78, 16, 3 },
{ 78, 22, 2 },
{ 78, 39, 2 },
{ 78, 42, 1 },
{ 78, 58, 1 },
{ 78, 59, 1 },
{ 78, 6, 3 },
{ 78, 64, 2 },
{ 78, 67, 1 },
{ 78, 76, 1 },
{ 78, 79, 1 },
{ 79, 16, 3 },
{ 79, 22, 2 },
{ 79, 3, 2 },
{ 79, 33, 3 },
{ 79, 37, 3 },
{ 79, 39, 2 },
{ 79, 48, 3 },
{ 79, 50, 4 },
{ 79, 55, 1 },
{ 79, 64, 2 },
{ 79, 70, 2 },
{ 79, 73, 1 },
{ 79, 7, 4 },
{ 79, 74, 21 },
{ 79, 78, 1 },
{ 79, 81, 3 },
{ 80, 16, 2 },
{ 80, 22, 2 },
{ 80, 39, 2 },
{ 80, 55, 1 },
{ 80, 64, 2 },
{ 81, 10, 2 },
{ 81, 16, 1 },
{ 81, 18, 2 },
{ 81, 20, 1 },
{ 81, 22, 2 },
{ 81, 24, 2 },
{ 81, 27, 2 },
{ 81, 28, 2 },
{ 81, 29, 2 },
{ 81, 33, 1 },
{ 81, 37, 1 },
{ 81, 39, 2 },
{ 81, 40, 1 },
{ 81, 42, 3 },
{ 81, 43, 2 },
{ 81, 48, 1 },
{ 81, 49, 1 },
{ 81, 50, 2 },
{ 81, 51, 2 },
{ 81, 5, 2 },
{ 81, 52, 1 },
{ 81, 53, 2 },
{ 81, 56, 1 },
{ 81, 57, 1 },
{ 81, 60, 2 },
{ 81, 61, 2 },
{ 81, 63, 2 },
{ 81, 64, 2 },
{ 81, 69, 2 },
{ 81, 7, 2 },
{ 81, 76, 3 },
{ 81, 79, 3 },
{ 81, 8, 2 },
{ 81, 82, 2 },
{ 81, 83, 2 },
{ 81, 86, 2 },
{ 82, 1, 1 },
{ 82, 12, 1 },
{ 82, 16, 2 },
{ 82, 18, 1 },
{ 82, 2, 1 },
{ 82, 21, 1 },
{ 82, 22, 3 },
{ 82, 26, 1 },
{ 82, 29, 1 },
{ 82, 32, 1 },
{ 82, 33, 2 },
{ 82, 37, 2 },
{ 82, 39, 3 },
{ 82, 41, 1 },
{ 82, 43, 1 },
{ 82, 45, 1 },
{ 82, 48, 2 },
{ 82, 50, 1 },
{ 82, 5, 1 },
{ 82, 52, 2 },
{ 82, 54, 1 },
{ 82, 6, 2 },
{ 82, 63, 1 },
{ 82, 64, 3 },
{ 82, 65, 1 },
{ 82, 68, 1 },
{ 82, 69, 1 },
{ 82, 7, 1 },
{ 82, 72, 1 },
{ 82, 75, 2 },
{ 82, 8, 1 },
{ 82, 81, 2 },
{ 83, 16, 1 },
{ 83, 20, 1 },
{ 83, 22, 2 },
{ 83, 33, 1 },
{ 83, 39, 2 },
{ 83, 64, 2 },
{ 83, 81, 1 },
{ 84, 16, 2 },
{ 84, 22, 1 },
{ 84, 3, 1 },
{ 84, 39, 1 },
{ 84, 42, 2 },
{ 84, 76, 1 },
{ 85, 16, 2 },
{ 85, 22, 3 },
{ 86, 28, 1 },
{ 86, 48, 2 },
{ 86, 5, 1 },
{ 86, 60, 1 },
{ 86, 6, 2 },
{ 86, 81, 2 },
{ 87, 41, 1 },
{ 87, 5, 1 },
{ 87, 63, 1 },
{ 87, 7, 1 },
{ 88, 16, 1 },
{ 88, 20, 1 },
{ 88, 39, 2 },
{ 88, 5, 2}
};
