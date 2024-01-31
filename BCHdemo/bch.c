#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef short int16_t;
typedef unsigned short uint16_t;
typedef unsigned char uint8_t;
typedef signed char int8_t;

#define len 200
#define len2  100

#define GOAL 63 // 需要生成的比特数
int x[GOAL * 2]; // bit位置

//int m = 6;
//uint16_t n = 63;
//int k = 7;
//int t = 15;
//int d = 31;

int m = 3;
uint16_t n = 7;
int k;
int t = 1;
int d;


int8_t p[6];
int8_t alpha_to[64];
int8_t index_of[64];
int8_t g[56];

int8_t row = 37;

uint8_t key_256[32];

uint8_t puf_binary_new[8 * 37];
uint8_t helper_data_new[7 * 37];
uint8_t key_per_row[37];


// 生成密钥
void genKey256() {
    // 获取密钥
    memset(puf_binary_new, 0, sizeof(puf_binary_new));

    int index = 0, i = 0;

    for (; i < 256; i++) {
        puf_binary_new[index++] = rand() % 255 + 1;
    }

    memset(key_256, 0, 32 * sizeof(uint8_t));
    memcpy(&key_256, &puf_binary_new, 32 * sizeof(uint8_t));
}

void genKeyPerRow(uint8_t *keys, uint8_t *result) {
    int index_result = 0, index_key = 0;
    uint8_t shift;
    for (int i = 0; i < 256; i++) {
        shift = 7 - (i % 8);
        if ((keys[index_key] >> shift) & 0x1) {
            result[index_result] = result[index_result] | 0x1;
        }
        if ((i + 1) % 7 == 0) {
            result[index_result] = result[index_result] << 1;
            index_result++;
        } else {
            result[index_result] = result[index_result] << 1;
        }
        if ((i + 1) % 8 == 0) {
            index_key++;
        }
    }
    result[index_result] = result[index_result] << 3;
}


void initialize_p() {
    int i;
    for (i = 1; i < m; i++)
        p[i] = 0;
    p[0] = p[m] = 1;
    if (m == 2) p[1] = 1;
    else if (m == 3) p[1] = 1;
    else if (m == 4) p[1] = 1;
    else if (m == 5) p[2] = 1;
    else if (m == 6) p[1] = 1;
    else if (m == 7) p[1] = 1;
    else if (m == 8) p[4] = p[5] = p[6] = 1;
    else if (m == 9) p[4] = 1;
    else if (m == 10) p[3] = 1;
    else if (m == 11) p[2] = 1;
    else if (m == 12) p[3] = p[4] = p[7] = 1;
    else if (m == 13) p[1] = p[3] = p[4] = 1;
    else if (m == 14) p[1] = p[11] = p[12] = 1;
    else if (m == 15) p[1] = 1;
    else if (m == 16) p[2] = p[3] = p[5] = 1;
    else if (m == 17) p[3] = 1;
    else if (m == 18) p[7] = 1;
    else if (m == 19) p[1] = p[5] = p[6] = 1;
    else if (m == 20) p[3] = 1;
    printf("p(x) = ");
    for (i = 0; i <= m; i++) {
        printf("%1d", p[i]);
    }
    printf("\n");

}

void generate_gf() {
    register int i, mask;

    mask = 1;
    alpha_to[m] = 0;
    for (i = 0; i < m; i++) {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (p[i] != 0)
            alpha_to[m] ^= mask;
        mask <<= 1;
    }
    index_of[alpha_to[m]] = m;
    mask >>= 1;
    for (i = m + 1; i < n; i++) {
        if (alpha_to[i - 1] >= mask)
            alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
        else
            alpha_to[i] = alpha_to[i - 1] << 1;
        index_of[alpha_to[i]] = i;
    }
    index_of[0] = -1;
}

void gen_poly() {
    register int8_t ii, jj, ll, kaux;
    register int8_t test, aux, nocycles, root, noterms, rdncy;
    int8_t cycle[(int) pow(2, m)][m], size[(int8_t) pow(2, m)], min[(int8_t) pow(2, m)], zeros[(int8_t) pow(2, m)];


    memset(cycle, 0, sizeof(cycle));
    memset(size, 0, sizeof(size));
    memset(min, 0, sizeof(min));
    memset(zeros, 0, sizeof(zeros));

    /* Generate cycle sets modulo n, n = 2**m - 1 */
    cycle[0][0] = 0;
    size[0] = 1;
    cycle[1][0] = 1;
    size[1] = 1;
    jj = 1;            /* cycle set index */
    do {
        /* Generate the jj-th cycle set */
        ii = 0;
        do {
            ii++;
            cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n;
            size[jj]++;
            aux = (cycle[jj][ii] * 2) % n;
        } while (aux != cycle[jj][0]);
        /* Next cycle set representative */
        ll = 0;
        do {
            ll++;
            test = 0;
            for (ii = 1; ((ii <= jj) && (!test)); ii++)
                /* Examine previous cycle sets */
                for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
                    if (ll == cycle[ii][kaux])
                        test = 1;
        } while ((test) && (ll < (n - 1)));
        if (!(test)) {
            jj++;    /* next cycle set index */
            cycle[jj][0] = ll;
            size[jj] = 1;
        }
    } while (ll < (n - 1));
    nocycles = jj;        /* number of cycle sets modulo n */

    d = 2 * t + 1;

    /* Search for roots 1, 2, ..., d-1 in cycle sets */
    kaux = 0;
    rdncy = 0;
    for (ii = 1; ii <= nocycles; ii++) {
        min[kaux] = 0;
        test = 0;
        for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
            for (root = 1; ((root < d) && (!test)); root++)
                if (root == cycle[ii][jj]) {
                    test = 1;
                    min[kaux] = ii;
                }
        if (min[kaux]) {
            rdncy += size[min[kaux]];
            kaux++;
        }
    }
    noterms = kaux;
    kaux = 1;
    for (ii = 0; ii < noterms; ii++)
        for (jj = 0; jj < size[min[ii]]; jj++) {
            zeros[kaux] = cycle[min[ii]][jj];
            kaux++;
        }

    k = n - rdncy;

    if (k < 0) {
        printf("Parameters invalid!\n");
        exit(0);
    }

    /* Compute the generator polynomial */
    g[0] = alpha_to[zeros[1]];
    g[1] = 1;        /* g(x) = (X + zeros[1]) initially */
    for (ii = 2; ii <= rdncy; ii++) {
        g[ii] = 1;
        for (jj = ii - 1; jj > 0; jj--)
            if (g[jj] != 0)
                g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
            else
                g[jj] = g[jj - 1];
        g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
    }


    printf("Generator polynomial:\ng(x) = ");
    for (ii = 0; ii <= rdncy; ii++) {
        printf("%d", g[ii]);
    }
    printf("\n");
}

void initBCH() {
    memset(p, 0, sizeof(p));
    memset(alpha_to, 0, sizeof(alpha_to));
    memset(index_of, 0, sizeof(index_of));
    memset(g, 0, sizeof(g));

    initialize_p();
    generate_gf();
    gen_poly();
}

void convertUint8ToBinArray(uint8_t input, int8_t *result) {
    for (int8_t i = 0; i < 8; i++) {
        if ((input >> (7 - i) & 0x1) == 1) {
            result[i] = 1;
        }
    }
}

void convertUint8ArrayToBinArray(uint8_t *input, int8_t *result, int length) {
    for (int8_t j = 0; j < length; j++) {
        convertUint8ToBinArray(input[j], &result[j * 8]);
    }
}

void convertBinArrayToUint8Array(int8_t *binary, uint8_t *result, int length) {
    int index = 0;
    int ro = length % 8 ? length + 8 - (length % 8) : length;
    for (int i = 0; i < ro; i++) {
        if (i < length) {
            if (binary[i] == 1) {
                result[index] = result[index] | 0x1;
            }
        }
        if ((i + 1) % 8 == 0) {
            index++;
        } else {
            result[index] = result[index] << 1;
        }
    }
}

void combineTwoBinArrayToUint8Array(int8_t *bin1, int length_bin1, int8_t *bin2, int length_bin2,
                                    uint8_t *result) {
    int8_t temp[length_bin1 + length_bin2];
    memcpy(temp, bin1, length_bin1);
    memcpy(&temp[length_bin1], bin2, length_bin2);
    convertBinArrayToUint8Array(temp, result, length_bin1 + length_bin2);
}

void encodeBch(uint8_t *in, uint8_t *result)
/*
 * Compute redundacy bb[], the coefficients of b(x). The redundancy
 * polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
 * by the generator polynomial g(x).
 */
{
    register int8_t i, j;
    register int8_t feedback;

    int8_t input[8];
    memset(input, 0, sizeof(input));
    convertUint8ArrayToBinArray(in, input, 1);

    int8_t bb[n - k];
    for (i = 0; i < n - k; i++)
        bb[i] = 0;
    for (i = k - 1; i >= 0; i--) {
        feedback = input[i] ^ bb[n - k - 1];
        if (feedback != 0) {
            for (j = n - k - 1; j > 0; j--)
                if (g[j] != 0)
                    bb[j] = bb[j - 1] ^ feedback;
                else
                    bb[j] = bb[j - 1];
            bb[0] = g[0] && feedback;
        } else {
            for (j = n - k - 1; j > 0; j--)
                bb[j] = bb[j - 1];
            bb[0] = 0;
        }
    }

    combineTwoBinArrayToUint8Array(bb, n - k, input, k, result);
}


void decodeBch(uint8_t *in, uint8_t *result) {
    int8_t i, j, u, q, t2, count = 0, syn_error = 0;
//    int8_t elp[(int) pow(2, m) + 2][(int) pow(2, m)], d[(int) pow(2, m) + 2], l[(int) pow(2, m) + 2], u_lu[
//        (int) pow(2, m) + 2], s[(int) pow(2, m) + 1];
    int8_t root[200], loc[200], err[(int) pow(2, m)], reg[201];

    int elp_m = (int) pow(2, m) + 2, elp_n = (int) pow(2, m);
    int8_t *elp = (int8_t *) malloc(sizeof(int8_t) * elp_m * elp_n),
            *d = (int8_t *) malloc(sizeof(int8_t) * ((int) pow(2, m) + 2)),
            *l = (int8_t *) malloc(sizeof(int8_t) * ((int) pow(2, m) + 2)),
            *u_lu = (int8_t *) malloc(sizeof(int8_t) * ((int) pow(2, m) + 2)),
            *s = (int8_t *) malloc(sizeof(int8_t) * ((int) pow(2, m) + 1));

    // memset(elp, 0, sizeof(elp));
    // memset(d, 0, sizeof(d));
    // memset(l, 0, sizeof(l));
    // memset(u_lu, 0, sizeof(u_lu));
    // memset(s, 0, sizeof(s));
    // memset(root, 0, sizeof(root));
    // memset(loc, 0, sizeof(loc));
    // memset(err, 0, sizeof(err));
    // memset(reg, 0, sizeof(reg));
    // memset(input, 0, sizeof(input));

    int8_t input[64];
    memset(input, 0, sizeof(input));
    convertUint8ArrayToBinArray(in, input, 8);

    // decode_bch_old(input, result);
    t2 = 2 * t;

    /* first form the syndromes */
    for (i = 1; i <= t2; i++) {
        s[i] = 0;
        for (j = 0; j < n; j++)
            if (input[j] != 0)
                s[i] ^= alpha_to[(i * j) % n];
        if (s[i] != 0)
            syn_error = 1; /* set error flag if non-zero syndrome */
        /*
         * Note:    If the code is used only for ERROR DETECTION, then
         *          exit program here indicating the presence of errors.
         */
        /* convert syndrome from polynomial form to index form  */
        s[i] = index_of[s[i]];
    }

    if (syn_error) {    /* if there are errors, try to correct them */
        /*
         * Compute the error location polynomial via the Berlekamp
         * iterative algorithm. Following the terminology of Lin and
         * Costello's book :   d[u] is the 'mu'th discrepancy, where
         * u='mu'+1 and 'mu' (the Greek letter!) is the step number
         * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
         * the elp at that step, and u_l[u] is the difference between
         * the step number and the degree of the elp.
         */
        /* initialise table entries */
        d[0] = 0;            /* index form */
        d[1] = s[1];        /* index form */
        *elp = 0;
        *(&elp + (int8_t) elp_n) = (int8_t *) 1;
//        elp[0][0] = 0;        /* index form */
//        elp[1][0] = 1;        /* polynomial form */
        for (i = 1; i < t2; i++) {
            *(&elp + i) = (int8_t *) -1;
            *(&elp + elp_n + i) = 0;
//            elp[0][i] = -1;    /* index form */
//            elp[1][i] = 0;    /* polynomial form */
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        do {
            u++;
            if (d[u] == -1) {
                l[u + 1] = l[u];
                for (i = 0; i <= l[u]; i++) {
                    *(&elp + (u + 1) * elp_n + i) = *(&elp + u * elp_n + i);
                    *(&elp + u * elp_n + i) = (int8_t *) (u * elp_n + i);
//                    elp[u + 1][i] = elp[u][i];
//                    elp[u][i] = index_of[elp[u][i]];
                }
            } else
                /*
                 * search for words with greatest u_lu[q] for
                 * which d[q]!=0
                 */
            {
                q = u - 1;
                while ((d[q] == -1) && (q > 0))
                    q--;
                /* have found first non-zero d[q]  */
                if (q > 0) {
                    j = q;
                    do {
                        j--;
                        if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
                            q = j;
                    } while (j > 0);
                }

                /*
                 * have now found q such that d[u]!=0 and
                 * u_lu[q] is maximum
                 */
                /* store degree of new elp polynomial */
                if (l[u] > l[q] + u - q)
                    l[u + 1] = l[u];
                else
                    l[u + 1] = l[q] + u - q;

                /* form new elp(x) */
                for (i = 0; i < t2; i++)
                    *(&elp + (u + 1) * elp_n + i) = 0;
//                    elp[u + 1][i] = 0;
                for (i = 0; i <= l[q]; i++)
                    if (*(&elp + q * elp_n + i) != -1)
                        *(&elp + (u + 1) * elp_n + i + u - q) = (int8_t *) alpha_to[
                                (d[u] + n - d[q] + (int8_t) (*(&elp + q * elp_n + i))) % n];
//                for (i = 0; i <= l[q]; i++)
//                    if (elp[q][i] != -1)
//                        elp[u + 1][i + u - q] =
//                            alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
                for (i = 0; i <= l[u]; i++) {
                    int8_t a = (int8_t) *(&elp + (u + 1) * elp_n + i);
                    int8_t b = (int8_t) *(&elp + u * elp_n + i);
                    *(&elp + (u + 1) * elp_n + i) = (int8_t *) (a ^ b);
                    *(&elp + u * elp_n + i) = (int8_t *) (u * elp_n + i);
//                    elp[u + 1][i] ^= elp[u][i];
//                    elp[u][i] = index_of[elp[u][i]];
                }
            }
            u_lu[u + 1] = u - l[u + 1];

            /* form (u+1)th discrepancy */
            if (u < t2) {
                /* no discrepancy computed on last iteration */
                if (s[u + 1] != -1)
                    d[u + 1] = alpha_to[s[u + 1]];
                else
                    d[u + 1] = 0;
                for (i = 1; i <= l[u + 1]; i++)
                    if ((s[u + 1 - i] != -1) && (*(&elp + (u + 1) * elp_n + 1 + i) != 0))
                        d[u + 1] ^= alpha_to[(s[u + 1 - i]
                                              + (u + 1) * elp_n + i) % n];
//                    if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
//                        d[u + 1] ^= alpha_to[(s[u + 1 - i]
//                                              + index_of[elp[u + 1][i]]) % n];
                /* put d[u+1] into index form */
                d[u + 1] = index_of[d[u + 1]];
            }
        } while ((u < t2) && (l[u + 1] <= t));

        u++;
        if (l[u] <= t) {/* Can correct errors */
            /* put elp into index form */
            for (i = 0; i <= l[u]; i++)
                *(&elp + u * elp_n + i) = (int8_t *) (u * elp_n + i);
//                elp[u][i] = index_of[elp[u][i]];

            /* Chien search: find roots of the error location polynomial */
            for (i = 1; i <= l[u]; i++)
                reg[i] = (int8_t) *(&elp + u * elp_n + i);
//                reg[i] = elp[u][i];
            count = 0;
            for (i = 1; i <= n; i++) {
                q = 1;
                for (j = 1; j <= l[u]; j++)
                    if (reg[j] != -1) {
                        reg[j] = (reg[j] + j) % n;
                        q ^= alpha_to[reg[j]];
                    }
                if (!q) {
                    /* store root and error
                    		 * location number indices */
                    root[count] = i;
                    loc[count] = n - i;
                    count++;
                }
            }
            if (count == l[u])
                /* no. roots = degree of elp hence <= t errors */
                for (i = 0; i < l[u]; i++)
                    input[loc[i]] ^= 1;
            // else    /* elp has degree >t hence cannot solve */
            //     Serial.println("Incomplete decoding: errors detected\n");
        }
    }

    convertBinArrayToUint8Array(&input[n - k], result, k);
//    memcpy(result, &input[n - k], sizeof(int8_t) * k);
}


void FuncOutputBin(uint8_t value) {
    char string[32] = {0};
    itoa(value, string, 2);
    printf("%s  ", string);
}

void main() {
    initBCH();

    printf("\nm = %d, n = %d t = %d k = %d ", m, n, t, k);

//
//    genKey256();
//
//    genKeyPerRow(key_256, key_per_row);

    // 7
    uint8_t key[1] = {
            0b01111100};//, 0b00001111, 0b10000110, 0b11000000, 0b01111000, 0b10010110, 0b01000011, 0b11111110};
    printf("\n\nkey: ");
    for (int i = 0; i < 1; i++) FuncOutputBin(key[i]);


    // 3
    uint8_t key_per = 0b10100000;
    printf("\n\nkey_per: ");
    FuncOutputBin(key_per);
    uint8_t encode[1];

    uint8_t xor[1];

    uint8_t helper_data[1];

//    uint8_t keyx1[8] = {0b00010010, 0b01001010,0b10000100,0b00101011,
//                        0b10010010,0b00110010,0b11000010,0b00110100
//                       };


    uint8_t keyx2[1] = {
            0b01101100};//, 0b00011111, 0b10000110, 0b11000000, 0b01111000, 0b10010110, 0b01000011, 0b11111110};

    uint8_t xorx[1];

    uint8_t keyx[1];

//    memcpy(encode, 0, sizeof(encode));

    encodeBch(&key_per, &encode[0]);
    printf("\n\nkey_per_encode: ");
    for (int i = 0; i < 1; i++) FuncOutputBin(encode[i]);

    for (int j = 0; j < 1; j++) {
        xor[j] = encode[j] ^ key[j];
    }

    memcpy(&helper_data[0], &xor[0], 8 * sizeof(uint8_t));

    printf("\n\nhelper_data: ");
    for (int j = 0; j < 1; j++) {
        FuncOutputBin(helper_data[j]);
    }

    for (int j = 0; j < 1; j++) {
        xorx[j] = helper_data[j] ^ keyx2[j];
    }
    printf("\n\nkeyx = ");
    for (int j = 0; j < 1; j++) {
        FuncOutputBin(keyx2[j]);
    }

    printf("\n\nxorx: ");
    for (int j = 0; j < 1; j++) {
        FuncOutputBin(xorx[j]);
    }

    decodeBch(&xorx, &keyx);

    printf("\n\nkeyx: ");
    for (int j = 0; j < 1; j++) {
        FuncOutputBin(keyx[j]);
    }

//    uint8_t encoded_new[8*37];
//    uint8_t xor_enroll_new[8*37];
//
//    int row = 37;
//    int n = n;
//    int k = k;

    // ASSERT ENCODED
//    for (int i = 0; i < row; i++) {
//        encodeBch(&key_per_row[i], &encoded_new[i * 8]);
//    }


    // ASSERT XOR ENROLL
//    for (int i = 0; i < row; i++) {
//        for (int j = 0; j < 8; j++) {
//            xor_enroll_new[i*8 + j] = encoded_new[i*8 + j] ^ puf_binary_new[i*8 + j];
//        }
//    }


    // ASSERT HELPER DATA
//    for (int i = 0; i < 37; i++) {
//        memcpy(&helper_data_new[i * 7], &xor_enroll_new[i * 8], 7 * sizeof(uint8_t));
//    }
//
//    for (int i = 0; i < 37*7; i++) {
//        printf("%d",helper_data_new[i]);
//    }
}
