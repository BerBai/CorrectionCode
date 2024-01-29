/* reedsolomon.h
 *
 * Reed-Solomon编码器/解码器的全局定义
 * 作者：Phil Karn（karn@ka9q.ampr.org），1996年9月
 * 版权所有 1999 年 Phil Karn, KA9Q
 *
 * Wireshark - 网络流量分析器
 * 作者：Gerald Combs <gerald@wireshark.org>
 * 版权所有 1998 年 Gerald Combs
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifdef __cplusplus
extern "C" {
#endif

/* 设置以下之一以启用编码器/解码器的调试和错误检查，但会影响速度 */
/* #undef DEBUG 1*/
/* #undef DEBUG 2*/
#undef DEBUG

/* 若要选择CCSDS标准（255,223）码，请定义CCSDS。这意味着MM、KK、B0和PRIM的标准值。*/
/* #undef CCSDS 1*/
#undef CCSDS
#ifndef CCSDS

/* 否则，保持CCSDS未定义，并设置以下参数：
 *
 * 将MM设置为每个代码符号的位数。Reed-Solomon分组大小将是NN = 2**M - 1个符号。支持的值在rs.c中定义。
 */
#define MM 8 /* 符号大小（位）*/

/*
 * 将KK设置为每个块中的数据符号数，必须小于块大小。代码将能够纠正最多NN-KK个擦除或（NN-KK）/2个错误，
 * 或者二者的组合，每个错误计为两个擦除。
 */
#define KK 207 /* 每个块中的数据符号数*/

/* 将B0设置为生成多项式的第一个根，以α形式表示，并将PRIM设置为用于生成生成多项式的根的α的幂。然后生成多项式将是
 * @**PRIM*B0, @**PRIM*(B0+1), @**PRIM*(B0+2)...@**PRIM*(B0+NN-KK)
 * “@”表示小写的α。
 */
#define B0 1 /* 生成多项式的第一个根，α形式*/
#define PRIM 1 /* 用于生成多项式的根的α的幂*/
#define STANDARD_ORDER

/* 如果您想选择自己的域生成多项式，您必须在rs.c中编辑它。*/

#else /* CCSDS */
/* 不要更改这些，它们是CCSDS标准 */
#define MM 8
#define KK 223
#define B0 112
#define PRIM 11
#endif

#define    NN ((1 << MM) - 1)

#if MM <= 8
typedef unsigned char dtype;
#else
typedef unsigned int dtype;
#endif

/* Reed-Solomon 编码
 * data[] 是输入块，纠错符号放在bb[]中
 * bb[] 可能位于数据块的末尾，例如对于（255，223）：
 *	encode_rs(&data[0],&data[223]);
 */
int encode_rs(dtype data[], dtype bb[]);

/* Reed-Solomon 擦除和错误诊断解码
 * 接收的块放入data[]中，零起始的擦除位置列表（如果有）放在eras_pos[]中，计数放在no_eras中。
 *
 * 如果可能，解码器会在原地纠正符号，并返回纠正的符号数。如果码字非法或不可纠正，则数据数组保持不变，返回-1。
 */
int eras_dec_rs(dtype data[], int eras_pos[], int no_eras);

#ifdef __cplusplus
}
#endif