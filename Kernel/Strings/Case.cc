#include<cstdio>
#include<cstring>
#include<cctype>
#include<cstdlib>

/*!
 * \details The original string is overwritten.
 * \param[in] str A C-string (char*), that will be overwritten.
 * \return str
 */
char* strtolower(char* str)
{
    for (int i = 0; i < strlen(str); i++)
        str[i] = tolower(str[i]);
    return str;
}

/*!
 * \details The original string is overwritten.
 * \param[in] str A C-string (char*), that will be overwritten.
 * \return str
 */
char* strtoupper(char* str)
{
    for (int i = 0; i < strlen(str); i++)
        str[i] = toupper(str[i]);
    return str;
}

/*!
 * \details The original strings are not overwritten.
 * \param[in] s1 A C-string (char*)
 * \param[in] s2 A C-string (char*)
 * \return As with strcmp, returns 0 if the strings are identical up to case, 
 * and otherwise 1.
 */
int strcicmp(const char* s1, const char* s2)
{
    char* t1 = (char*) calloc(strlen(s1) + 1, sizeof(char));
    char* t2 = (char*) calloc(strlen(s2) + 1, sizeof(char));
    strncpy(t1, s1, strlen(s1));
    strncpy(t2, s2, strlen(s2));
    strtolower(t1);
    strtolower(t2);
    int out = strcmp(t1, t2);
    if (out != 0)
        out = 1;
    
    free(t1);
    free(t2);
    return out;
}
