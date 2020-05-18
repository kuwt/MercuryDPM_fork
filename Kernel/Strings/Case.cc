//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
