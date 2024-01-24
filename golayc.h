/*
golay.h: 实现二进制Golay码[24,12,8]

作者: Axel Persinger
许可证: MIT许可证
*/

#ifndef _GOLAY_H
#define _GOLAY_H

/*
包含语句

stdio - 标准I/O
stdlib - 标准库
*/
#include <stdio.h>
#include <stdlib.h>


/*
预处理器定义

POLY - Golay的特征多项式（您可以选择0xAE3或0xC75，都可以）
*/
#define POLY 0xAE3


/*
类型定义

codeword - 存储单个码字的所有信息
*/
typedef struct codeword {
    unsigned short parity: 1;
    unsigned short check: 12;
    unsigned short information: 12;
} st_codeword;


/*
函数声明

add - 简单的CFFI示例
calculate_check - 计算并填写码字的检验位
calculate_parity - 计算并填写码字的奇偶校验位
*/
int add(int i, int j);
void calculate_check(st_codeword* code);
void calculate_parity(st_codeword* code);


#endif