#include <stdio.h>
#include <string.h>
#include <stdbool.h>

int main() {
    char str1[] = "eq"; // "abcdef";
    char str2[] = "eq";
    int result;
    bool is_state_vector, is_ecliptic = true;

    
    result = memcmp(str1, str2, 2); // Compare the first 3 bytes
    printf("Comparison result: %d\n", result);
    if( !memcmp( str1, "eq", 2))
         is_ecliptic = false;
    printf("is_ecliptic=%d\n", is_ecliptic);
    if ( !is_ecliptic)
    {
        printf("Calling equatorial_to_ecliptic\n");
    }
    result = memcmp(str1, str2, 4); // Compare the first 4 bytes
    printf("Comparison result: %d\n", result);

    return 0;
}
