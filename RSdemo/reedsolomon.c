/* reedsolomon.c
 *
 * Reed-Solomon编码和解码，由Phil Karn（karn@ka9q.ampr.org）在1996年9月创建
 * 版权所有 1999 年 Phil Karn, KA9Q
 * 单独的CCSDS版本在1998年12月创建，于1999年5月合并到此版本
 *
 * 该文件源自我的通用RS编码器/解码器，该编码器/解码器又基于Robert Morelos-Zaragoza（robert@spectra.eng.hawaii.edu）和Hari Thirumoorthy（harit@spectra.eng.hawaii.edu）于1995年8月创作的“new_rs_erasures.c”程序。
 *
 * Wireshark - 网络流量分析器
 * 作者：Gerald Combs <gerald@wireshark.org>
 * 版权所有 1998 年 Gerald Combs
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "reedsolomon.h"
//#include <wsutil/ws_printf.h> /* ws_debug_printf */

#ifdef CCSDS
/* CCSDS field generator polynomial: 1+x+x^2+x^7+x^8 */
int Pp[MM+1] = { 1, 1, 1, 0, 0, 0, 0, 1, 1 };

#else /* not CCSDS */
/* MM, KK, B0, PRIM are user-defined in rs.h */

/* Primitive polynomials - see Lin & Costello, Appendix A,
 * and  Lee & Messerschmitt, p. 453.
 */
#if(MM == 2)/* Admittedly silly */
int Pp[MM+1] = { 1, 1, 1 };

#elif(MM == 3)
/* 1 + x + x^3 */
int Pp[MM+1] = { 1, 1, 0, 1 };

#elif(MM == 4)
/* 1 + x + x^4 */
int Pp[MM+1] = { 1, 1, 0, 0, 1 };

#elif(MM == 5)
/* 1 + x^2 + x^5 */
int Pp[MM+1] = { 1, 0, 1, 0, 0, 1 };

#elif(MM == 6)
/* 1 + x + x^6 */
int Pp[MM+1] = { 1, 1, 0, 0, 0, 0, 1 };

#elif(MM == 7)
/* 1 + x^3 + x^7 */
int Pp[MM+1] = { 1, 0, 0, 1, 0, 0, 0, 1 };

#elif(MM == 8)
/* 1+x^2+x^3+x^4+x^8 */
int Pp[MM + 1] = {1, 0, 1, 1, 1, 0, 0, 0, 1};

#elif(MM == 9)
/* 1+x^4+x^9 */
int Pp[MM+1] = { 1, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

#elif(MM == 10)
/* 1+x^3+x^10 */
int Pp[MM+1] = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

#elif(MM == 11)
/* 1+x^2+x^11 */
int Pp[MM+1] = { 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

#elif(MM == 12)
/* 1+x+x^4+x^6+x^12 */
int Pp[MM+1] = { 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 };

#elif(MM == 13)
/* 1+x+x^3+x^4+x^13 */
int Pp[MM+1] = { 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

#elif(MM == 14)
/* 1+x+x^6+x^10+x^14 */
int Pp[MM+1] = { 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };

#elif(MM == 15)
/* 1+x+x^15 */
int Pp[MM+1] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

#elif(MM == 16)
/* 1+x+x^3+x^12+x^16 */
int Pp[MM+1] = { 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 };

#else
#error "Either CCSDS must be defined, or MM must be set in range 2-16"
#endif

#endif

#ifdef STANDARD_ORDER /* first byte transmitted is index of x**(KK-1) in message poly*/
/* definitions used in the encode routine*/
#define MESSAGE(i) data[KK-(i)-1]
#define REMAINDER(i) bb[NN-KK-(i)-1]
/* definitions used in the decode routine*/
#define RECEIVED(i) data[NN-1-(i)]
#define ERAS_INDEX(i) (NN-1-eras_pos[i])
#define INDEX_TO_POS(i) (NN-1-(i))
#else /* first byte transmitted is index of x**0 in message polynomial*/
/* definitions used in the encode routine*/
#define MESSAGE(i) data[i]
#define REMAINDER(i) bb[i]
/* definitions used in the decode routine*/
#define RECEIVED(i) data[i]
#define ERAS_INDEX(i) eras_pos[i]
#define INDEX_TO_POS(i) i
#endif


/* 此处定义了用于存储代码中使用的有限域元素的类型。
 * 请确保此类型大于 char，如果使用的是 GF(256) 以上的值。
 *
 * 注意：无符号字符类型可用于 GF(256)，但在奔腾处理器上，int 似乎更快。
 */
typedef int gf;

/* index->polynomial form conversion table */
static gf Alpha_to[NN + 1];

/* Polynomial->index form conversion table */
static gf Index_of[NN + 1];

/* No legal value in index form represents zero, so
 * we need a special value for this purpose
 */
#define A0      (NN)

/* Generator polynomial g(x) in index form */
static gf Gg[NN - KK + 1];

static int RS_init; /* Initialization flag */

/* 计算 x % NN，其中 NN 是 2 ** MM - 1，
 * 无需缓慢的除法
 */
/* static inline gf*/
static gf modnn(int x) {
    while (x >= NN) {
        x -= NN;
        x = (x >> MM) + (x & NN);
    }
    return x;
}

#define min_(a, b)       ((a) < (b) ? (a) : (b))

#define CLEAR(a, n) {\
int ci;\
for(ci=(n)-1;ci >=0;ci--)\
(a)[ci] = 0;\
}

#define COPY(a, b, n) {\
int ci;\
for(ci=(n)-1;ci >=0;ci--)\
(a)[ci] = (b)[ci];\
}

#define COPYDOWN(a, b, n) {\
int ci;\
for(ci=(n)-1;ci >=0;ci--)\
(a)[ci] = (b)[ci];\
}

static void init_rs(void);

#ifdef CCSDS
/* 从常规alpha转换为Berlekamp的双基表示的转换查找表。仅在CCSDS版本中使用。
 * taltab[] -- 将常规转换为双基
 * tal1tab[] -- 将双基转换为常规
 *
 * 注意：实际的RS编码器/解码器使用常规基础。因此，在编码或解码之前，数据会从双基转换为常规基，然后再转换回来。
 */
static unsigned char taltab[NN+1],tal1tab[NN+1];

static unsigned char tal[] = { 0x8d, 0xef, 0xec, 0x86, 0xfa, 0x99, 0xaf, 0x7b };

/* 生成常规alpha表示和Berlekamp的双基表示之间的转换查找表
 * (@**7, @**6, ...@**0)
 *  and Berlekamp's dual basis representation
 * (l0, l1, ...l7)
 */
static void gen_ltab(void)
{
  int i,j,k;

  for(i=0;i<256;i++){/* 对于每个输入值 */
    taltab[i] = 0;
    for(j=0;j<8;j++)/* 对于矩阵的每一列 */
      for(k=0;k<8;k++){/* 对于矩阵的每一行 */
        if(i & (1<<k))
          taltab[i] ^= tal[7-k] & (1<<j);
      }
    tal1tab[taltab[i]] = i;
  }
}
#endif /* CCSDS */

#if PRIM != 1
static int Ldec;/* Decrement for aux location variable in Chien search */

static void gen_ldec(void)
{
    for(Ldec=1;(Ldec % PRIM) != 0;Ldec+= NN)
        ;
    Ldec /= PRIM;
}
#else
#define Ldec 1
#endif

/* 从Pp[0]..Pp[m]的不可约多项式p(X)生成GF(2**m)
 * 查找表：索引->多项式形式   alpha_to[] 包含 j=alpha**i;
 * 多项式形式 -> 索引形式   index_of[j=alpha**i] = i
其中alpha=2是GF(2**m)的原始元素
HARI的评论: (4/13/94) 可以如下使用alpha_to[]：
假设@表示通常称为"alpha"的原始元素，即原始多项式p(x)的根。那么在GF(2^m)中，对于任何 0 <= i <= 2^m-2，
@^i = a(0) + a(1) @ + a(2) @^2 + ... + a(m-1) @^(m-1)
其中二进制向量(a(0), a(1), a(2),..., a(m-1))是整数"alpha_to[i]"的表示，a(0)是LSB，a(m-1)是MSB。因此，例如@^5的多项式表示将由整数"alpha_to[5]"的二进制表示给出。
类似地，index_of[]可以使用如下：
如上所述，让@表示GF(2^m)的原始元素，即原始多项式p(x)的根。为了找到@（alpha）的幂，其多项式表示为
a(0) + a(1) @ + a(2) @^2 + ... + a(m-1) @^(m-1) ，我们考虑整数"i"，其二进制表示法是(a(0)作为LSB和a(m-1)作为MSB)，然后定位条目
"index_of[i]"。现在，@^index_of[i]是那个元素，其多项式表示为(a(0),a(1),a(2),...,a(m-1))。
注意：
元素alpha_to[2^m-1] = 0 总是表示"@^infinity" = 0 的表示为(0,0,0,...,0)。
类似地，元素index_of[0] = A0 总是表示具有多项式表示(0,0,...,0)的alpha的幂是“无穷大”。

*/

static void generate_gf(void) {
    register int i, mask;

    mask = 1;
    Alpha_to[MM] = 0;
    for (i = 0; i < MM; i++) {
        Alpha_to[i] = mask;
        Index_of[Alpha_to[i]] = i;
        /* If Pp[i] == 1 then, term @^i occurs in poly-repr of @^MM */
        if (Pp[i] != 0)
            Alpha_to[MM] ^= mask;     /* Bit-wise EXOR operation */
        mask <<= 1; /* single left-shift */
    }
    Index_of[Alpha_to[MM]] = MM;
    /*
     * 已获得@^MM的多项式表示。
     * @^(i+1)的多项式表示由@^i的多项式表示向左移动一个比特，
     * 并考虑@^i被移位时可能出现的任何@^MM项。
     */
    mask >>= 1;
    for (i = MM + 1; i < NN; i++) {
        if (Alpha_to[i - 1] >= mask)
            Alpha_to[i] = Alpha_to[MM] ^ ((Alpha_to[i - 1] ^ mask) << 1);
        else
            Alpha_to[i] = Alpha_to[i - 1] << 1;
        Index_of[Alpha_to[i]] = i;
    }
    Index_of[0] = A0;
    Alpha_to[NN] = 0;
}

/* 从（X+@**(B0+i)）的乘积中获取TT错误更正，长度为NN=(2**MM -1)的Reed Solomon码的生成多项式。

 * 例子：

 * 如果B0 = 1，TT = 1。g(x)的次数 = 2*TT = 2。
 * g(x) = (x+@) (x+@**2)

 * 如果B0 = 0，TT = 2。g(x)的次数 = 2*TT = 4。
 * g(x) = (x+1) (x+@) (x+@**2) (x+@**3)
 **/
static void gen_poly(void) {
    register int i, j;

    Gg[0] = 1;
    for (i = 0; i < NN - KK; i++) {
        Gg[i + 1] = 1;
        /*
         * Below multiply (Gg[0]+Gg[1]*x + ... +Gg[i]x^i) by
         * (@**(B0+i)*PRIM + x)
         */
        for (j = i; j > 0; j--)
            if (Gg[j] != 0)
                Gg[j] = Gg[j - 1] ^ Alpha_to[modnn((Index_of[Gg[j]]) + (B0 + i) * PRIM)];
            else
                Gg[j] = Gg[j - 1];
        /* Gg[0] can never be zero */
        Gg[0] = Alpha_to[modnn(Index_of[Gg[0]] + (B0 + i) * PRIM)];
    }
    /* convert Gg[] to index form for quicker encoding */
    for (i = 0; i <= NN - KK; i++)
        Gg[i] = Index_of[Gg[i]];
}


/*
 * 对数据数组中的符号字符串（i=0..(k-1)）进行系统编码，产生 NN-KK 个校验符号，
 * 存储在 bb[0]..bb[NN-KK-1]，其中 data[] 是输入，bb[] 是以多项式形式输出。
 * 编码是通过使用由 Gg[] 元素指定的适当连接的反馈移位寄存器来完成的，Gg[] 是在上面生成的。
 * 编码字为 c(X) = data(X)*X**(NN-KK) + b(X)
 */

int encode_rs(dtype data[KK], dtype bb[NN - KK]) {
    register int i, j;
    gf feedback;

#if DEBUG >= 1 && MM != 8
    /* Check for illegal input values */
  for(i=0;i<KK;i++)
    if(MESSAGE(i) > NN)
      return -1;
#endif

    if (!RS_init)
        init_rs();

    CLEAR(bb, NN - KK);

#ifdef CCSDS
    /* Convert to conventional basis */
  for(i=0;i<KK;i++)
    MESSAGE(i) = tal1tab[MESSAGE(i)];
#endif

    for (i = KK - 1; i >= 0; i--) {
        feedback = Index_of[MESSAGE(i) ^ REMAINDER(NN - KK - 1)];
        if (feedback != A0) {       /* feedback term is non-zero */
            for (j = NN - KK - 1; j > 0; j--)
                if (Gg[j] != A0)
                    REMAINDER(j) = REMAINDER(j - 1) ^ Alpha_to[modnn(Gg[j] + feedback)];
                else
                    REMAINDER(j) = REMAINDER(j - 1);
            REMAINDER(0) = Alpha_to[modnn(Gg[0] + feedback)];
        } else {    /* feedback term is zero. encoder becomes a
                 * single-byte shifter */
            for (j = NN - KK - 1; j > 0; j--)
                REMAINDER(j) = REMAINDER(j - 1);
            REMAINDER(0) = 0;
        }
    }
#ifdef CCSDS
    /* Convert to l-basis */
  for(i=0;i<NN;i++)
    MESSAGE(i) = taltab[MESSAGE(i)];
#endif

    return 0;
}

/* 对RS码进行错误和擦除解码。如果解码成功，将编码字写入data[]本身。否则，data[]保持不变。
 * 返回更正的符号数量，如果编码字不合法或不可更正则返回-1。如果eras_pos非空，则将检测到的错误位置写回。
 * 注意！该数组的长度必须至少为NN-KK个元素。
 * 首先，由调用程序声明no_eras个擦除。
 * 然后，可更正的最大错误数量为t_after_eras = floor((NN-KK-no_eras)/2)。
 * 如果信道错误数量不超过"t_after_eras"，则传输的编码字将被恢复。
 * 算法的详细信息可以在R. Blahut的"Error-Correcting Codes理论"中找到。
 * 警告：eras_pos[]数组不得包含重复条目；解码器将失败。
 * 解码器*可以*检查此条件，但这会增加每次解码操作的额外时间。
 * */

int eras_dec_rs(dtype data[NN], int eras_pos[NN - KK], int no_eras) {
    int deg_lambda, el, deg_omega;
    int i, j, r, k;
    gf u, q, tmp, num1, num2, den, discr_r;
    gf lambda[NN - KK + 1], s[NN - KK + 1];   /* Err+Eras Locator poly
                                         * and syndrome poly */
    gf b[NN - KK + 1], t[NN - KK + 1], omega[NN - KK + 1];
    gf root[NN - KK], reg[NN - KK + 1], loc[NN - KK];
    int syn_error, count;

    if (!RS_init)
        init_rs();

#ifdef CCSDS
    /* Convert to conventional basis */
  for(i=0;i<NN;i++)
    RECEIVED(i) = tal1tab[RECEIVED(i)];
#endif

#if DEBUG >= 1 && MM != 8
    /* Check for illegal input values */
  for(i=0;i<NN;i++)
    if(RECEIVED(i) > NN)
      return -1;
#endif
    /* form the syndromes; i.e., evaluate data(x) at roots of g(x)
     * namely @**(B0+i)*PRIM, i = 0, ... ,(NN-KK-1)
     */
    for (i = 1; i <= NN - KK; i++) {
        s[i] = RECEIVED(0);
    }
    for (j = 1; j < NN; j++) {
        if (RECEIVED(j) == 0)
            continue;
        tmp = Index_of[RECEIVED(j)];

        /*  s[i] ^= Alpha_to[modnn(tmp + (B0+i-1)*j)]; */
        for (i = 1; i <= NN - KK; i++)
            s[i] ^= Alpha_to[modnn(tmp + (B0 + i - 1) * PRIM * j)];
    }
    /* Convert syndromes to index form, checking for nonzero condition */
    syn_error = 0;
    for (i = 1; i <= NN - KK; i++) {
        syn_error |= s[i];
        /*ws_debug_printf("syndrome %d = %x\n",i,s[i]);*/
        s[i] = Index_of[s[i]];
    }

    if (!syn_error) {
        /* if syndrome is zero, data[] is a codeword and there are no
         * errors to correct. So return data[] unmodified
         */
        count = 0;
        goto finish;
    }
    CLEAR(&lambda[1], NN - KK);
    lambda[0] = 1;

    if (no_eras > 0) {
        /* Init lambda to be the erasure locator polynomial */
        lambda[1] = Alpha_to[modnn(PRIM * ERAS_INDEX(0))];
        for (i = 1; i < no_eras; i++) {
            u = modnn(PRIM * ERAS_INDEX(i));
            for (j = i + 1; j > 0; j--) {
                tmp = Index_of[lambda[j - 1]];
                if (tmp != A0)
                    lambda[j] ^= Alpha_to[modnn(u + tmp)];
            }
        }
#if DEBUG >= 1
        /* Test code that verifies the erasure locator polynomial just constructed
       Needed only for decoder debugging. */

    /* find roots of the erasure location polynomial */
    for(i=1;i<=no_eras;i++)
      reg[i] = Index_of[lambda[i]];
    count = 0;
    for (i = 1,k=NN-Ldec; i <= NN; i++,k = modnn(NN+k-Ldec)) {
      q = 1;
      for (j = 1; j <= no_eras; j++)
        if (reg[j] != A0) {
          reg[j] = modnn(reg[j] + j);
          q ^= Alpha_to[reg[j]];
        }
      if (q != 0)
        continue;
      /* store root and error location number indices */
      root[count] = i;
      loc[count] = k;
      count++;
    }
    if (count != no_eras) {
      ws_debug_printf("\n lambda(x) is WRONG\n");
      count = -1;
      goto finish;
    }
#if DEBUG >= 2
    ws_debug_printf("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
    for (i = 0; i < count; i++)
      ws_debug_printf("%d ", loc[i]);
    ws_debug_printf("\n");
#endif
#endif
    }
    for (i = 0; i < NN - KK + 1; i++)
        b[i] = Index_of[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = no_eras;
    el = no_eras;
    while (++r <= NN - KK) {        /* r is the step number */
        /* Compute discrepancy at the r-th step in poly-form */
        discr_r = 0;
        for (i = 0; i < r; i++) {
            if ((lambda[i] != 0) && (s[r - i] != A0)) {
                discr_r ^= Alpha_to[modnn(Index_of[lambda[i]] + s[r - i])];
            }
        }
        discr_r = Index_of[discr_r];        /* Index form */
        if (discr_r == A0) {
            /* 2 lines below: B(x) <-- x*B(x) */
            COPYDOWN(&b[1], b, NN - KK);
            b[0] = A0;
        } else {
            /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
            t[0] = lambda[0];
            for (i = 0; i < NN - KK; i++) {
                if (b[i] != A0)
                    t[i + 1] = lambda[i + 1] ^ Alpha_to[modnn(discr_r + b[i])];
                else
                    t[i + 1] = lambda[i + 1];
            }
            if (2 * el <= r + no_eras - 1) {
                el = r + no_eras - el;
                /*
                 * 2 lines below: B(x) <-- inv(discr_r) *
                 * lambda(x)
                 */
                for (i = 0; i <= NN - KK; i++)
                    b[i] = (lambda[i] == 0) ? A0 : modnn(Index_of[lambda[i]] - discr_r + NN);
            } else {
                /* 2 lines below: B(x) <-- x*B(x) */
                COPYDOWN(&b[1], b, NN - KK);
                b[0] = A0;
            }
            COPY(lambda, t, NN - KK + 1);
        }
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for (i = 0; i < NN - KK + 1; i++) {
        lambda[i] = Index_of[lambda[i]];
        if (lambda[i] != A0)
            deg_lambda = i;
    }
    /*
     * Find roots of the error+erasure locator polynomial by Chien
     * Search
     */
    COPY(&reg[1], &lambda[1], NN - KK);
    count = 0;            /* Number of roots of lambda(x) */
    for (i = 1, k = NN - Ldec; i <= NN; i++, k = modnn(NN + k - Ldec)) {
        q = 1;
        for (j = deg_lambda; j > 0; j--) {
            if (reg[j] != A0) {
                reg[j] = modnn(reg[j] + j);
                q ^= Alpha_to[reg[j]];
            }
        }
        if (q != 0)
            continue;
        /* store root (index-form) and error location number */
        root[count] = i;
        loc[count] = k;
        /* If we've already found max possible roots,
         * abort the search to save time
         */
        if (++count == deg_lambda)
            break;
    }
    if (deg_lambda != count) {
        /*
         * deg(lambda) unequal to number of roots => uncorrectable
         * error detected
         */
        count = -1;
        goto finish;
    }
    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**(NN-KK)). in index form. Also find deg(omega).
     */
    deg_omega = 0;
    for (i = 0; i < NN - KK; i++) {
        tmp = 0;
        j = (deg_lambda < i) ? deg_lambda : i;
        for (; j >= 0; j--) {
            if ((s[i + 1 - j] != A0) && (lambda[j] != A0))
                tmp ^= Alpha_to[modnn(s[i + 1 - j] + lambda[j])];
        }
        if (tmp != 0)
            deg_omega = i;
        omega[i] = Index_of[tmp];
    }
    omega[NN - KK] = A0;

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(B0-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    for (j = count - 1; j >= 0; j--) {
        num1 = 0;
        for (i = deg_omega; i >= 0; i--) {
            if (omega[i] != A0)
                num1 ^= Alpha_to[modnn(omega[i] + i * root[j])];
        }
        num2 = Alpha_to[modnn(root[j] * (B0 - 1) + NN)];
        den = 0;

        /* lambda[i+1] for i even is the formal derivative lambda_pr of lambda[i] */
        for (i = min_(deg_lambda, NN - KK - 1) & ~1; i >= 0; i -= 2) {
            if (lambda[i + 1] != A0)
                den ^= Alpha_to[modnn(lambda[i + 1] + i * root[j])];
        }
        if (den == 0) {
#if DEBUG >= 1
            ws_debug_printf("\n ERROR: denominator = 0\n");
#endif
            /* Convert to dual- basis */
            count = -1;
            goto finish;
        }
        /* Apply error to data */
        if (num1 != 0) {
            RECEIVED(loc[j]) ^= Alpha_to[modnn(Index_of[num1] + Index_of[num2] + NN - Index_of[den])];
        }
    }
    finish:
#ifdef CCSDS
    /* Convert to dual- basis */
  for(i=0;i<NN;i++)
    RECEIVED(i) = taltab[RECEIVED(i)];
#endif
    if (eras_pos != NULL) {
        for (i = 0; i < count; i++) {
            if (eras_pos != NULL)
                eras_pos[i] = INDEX_TO_POS(loc[i]);
        }
    }
    return count;
}

/* Encoder/decoder initialization - call this first! */
static void init_rs(void) {
    generate_gf();
    gen_poly();
#ifdef CCSDS
    gen_ltab();
#endif
#if PRIM != 1
    gen_ldec();
#endif
    RS_init = 1;
}

/*
 * Editor modelines  -  https://www.wireshark.org/tools/modelines.html
 *
 * Local Variables:
 * c-basic-offset: 2
 * tab-width: 8
 * indent-tabs-mode: nil
 * End:
 *
 * ex: set shiftwidth=2 tabstop=8 expandtab:
 * :indentSize=2:tabSize=8:noTabs=true:
 */


void main() {
    system("chcp 65001");

    printf("MM: %d, NN: %d, KK: %d\n", MM, NN, KK);

    clock_t start, end;
    start = clock();
    init_rs();
    dtype data[NN] = {1, 2, 3, 4, 5, 6, 7, 8};
    int en = encode_rs(&data[0], &data[KK]);
    printf("%d\n", en);
    printf("信息：");
    for (int i = 0; i < NN; i++) {
        printf("%x ", data[i]);
    }
    printf("\n");
    for (int i = 0; i < 16; i++)
        data[i + 10] = 1 + i;
    int de = eras_dec_rs(&data[0], &data[KK], 0);

    end = clock();

    printf("检测错误数：%d\n", de);
    printf("纠错：");
    for (int i = 0; i < NN; i++) {
        printf("%x ", data[i]);
    }
    printf("\n");
    printf("用时：%lfms\n", (double) (end - start) * 1000 / CLOCKS_PER_SEC);
}