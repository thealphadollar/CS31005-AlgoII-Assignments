/*
Name: Shivam Kumar Jha
Roll: 17CS30033
*/

#include <math.h>
#include <stdio.h>


struct Coordinate
{
	float x, y;
};

float calcDist(struct Coordinate a, struct Coordinate b)
{
	return sqrt(((a.x - b.x)*(a.x - b.x)) + ((a.y - b.y)*(a.y - b.y)));
}

int divideOut(struct Coordinate arr[], int start, int last, char dir){
    if (dir == 'X')
    {
        int pivot = arr[last].y;    
        int i = (start - 1);  
        struct Coordinate temp;
        for (int j = start; j <= last- 1; j++) 
        {     
            if (arr[j].y <= pivot) 
            { 
                i++;    
                temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            } 
        } 

        temp = arr[i+1];
        arr[i+1] = arr[last];
        arr[last] = temp;
    
        return (i + 1);
    }
}

int divideIn(struct Coordinate arr[], int start, int last, char dir) 
{ 
    if (dir == 'Y')
    {
        int pivot = arr[last].y;    
        int i = (start - 1);  
        struct Coordinate temp;
        for (int j = start; j <= last- 1; j++) 
        {     
            if (arr[j].y <= pivot) 
            { 
                i++;    
                temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            } 
        } 

        temp = arr[i+1];
        arr[i+1] = arr[last];
        arr[last] = temp;
    
        return (i + 1);
    }

    else
    {
        int pivot = arr[last].x;    
        int i = (start - 1);  
        struct Coordinate temp;
        for (int j = start; j <= last- 1; j++) 
        {
            if (arr[j].x <= pivot) 
            { 
                i++;    
                temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            } 
        } 

        temp = arr[i+1];
        arr[i+1] = arr[last];
        arr[last] = temp;
    
        return i + 1; 
    }
}

void verticalSort(struct Coordinate arr[], int begin, int last) 
{ 
    if (begin >= last) return;
    int parti = divideIn(arr, begin, last, 'Y');   
    verticalSort(arr, begin, parti - 1); 
    verticalSort(arr, parti + 1, last);
}

void horziSort(struct Coordinate arr[], int begin, int last) 
{ 
    if (begin >= last) return;      
    int parti = divideIn(arr, begin, last, 'X'); 
    horziSort(arr, begin, parti - 1); 
    horziSort(arr, parti + 1, last); 
}

float findCloestInVertical(struct Coordinate *C, int n, float min_distance)
{
	int j = 0, mid = n/2;
	struct Coordinate mid_point = C[mid], cur_band[n];

	for (int i = 0; i < n; i++)
	{
		if (fabs(C[i].x - mid_point.x) < min_distance)
		{
			cur_band[j] = C[i];
			j++;
		}
	}

	verticalSort(cur_band, 0, j-1);

	float min = min_distance;
	for (int i = 0; i < j; ++i) 
	{
        for (int k = i+1; ((k < j) && ((cur_band[k].y - cur_band[i].y) < min)); ++k) 
        {	
            if (calcDist(cur_band[i], cur_band[k]) < min) 
               	 min = calcDist(cur_band[i], cur_band[k]); 
		}
	}

	if (min > min_distance)
		return min_distance;
    return min;
}


float findClosestPair(struct Coordinate *C, int n, int m) // m has the total number of points
{
	if (n <= 3)
	{
		float cur_min = 200*m + 1; 
    	for (int i = 0; i < n; ++i)
    	{
    		for (int j = i + 1; j < n; ++j)
    		{
    			if (calcDist(C[i], C[j]) < cur_min)
    			{
    				cur_min = calcDist(C[i],C[j]);
    			}
    		}
    	}
    return cur_min; 
	}

	int mid = n/2;
	struct Coordinate middle_point = C[mid];

	float left_plane = findClosestPair(C, mid, m);
	float right_plane = findClosestPair(C + mid, n - mid, m);
	float cur_min_distance;

	if (left_plane <= right_plane)
	{
		cur_min_distance = left_plane;
	}
	else
	{
		cur_min_distance = right_plane;
	}

	return findCloestInVertical(C, n, cur_min_distance);
}

int main()
{
	int num_points; // Number of points
	printf("\nNumber of points: ");
	scanf("%d", &num_points);
	struct Coordinate plane[num_points];

	for (int i = 0; i < num_points; i++)
	{
		plane[i].x = (101 + rand()%(100*num_points));
		plane[i].y = (101 + rand()%(100*num_points));
	}

    horziSort(plane, 0, num_points-1);
	float min_dist = findClosestPair(plane, num_points, num_points);

	FILE *file_ptr = fopen("result.svg", "w");
	fprintf(file_ptr, "<svg xmlns=\"http://www.w3.org/2000/svg\"><rect width=\"%d\" height=\"%d\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />", 144*num_points, 144*num_points);
	for (int i = 0; i < num_points; i++)
	{
		fprintf(file_ptr, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"#000000\" stroke-width=\"0\" fill=\"#0066ff\" fill-opacity=\"0.4\" /><circle cx=\"%f\" cy=\"%f\" r=\"2\" stroke=\"black\" stroke-width=\"0\" fill=\"#000000\" /><text x = \"%f\" y = \"%f\"> %d </text>", plane[i].x, plane[i].y, min_dist/2, plane[i].x, plane[i].y, plane[i].x, plane[i].y, i + 1);
	}
	fprintf(file_ptr, "</svg>");
	fclose(file_ptr);

    printf("SVG image generated to result.svg\n");
	return 0;
}