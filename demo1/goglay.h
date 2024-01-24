//
// Created by iaciac on 2024/1/24.
//

#ifndef GOLAYCODE_GOGLAY_H
#define GOLAYCODE_GOGLAY_H

#include <stdio.h>

#define POLY  0xAE3  /* 或者使用其他多项式，0xC75 */

unsigned long golay(unsigned long cw)
/* 此函数计算[23,12] Golay码字。
   返回的长整型格式为
   [校验位(11),数据(12)]。 */
{
    int i;
    unsigned long c;
    cw &= 0xfffl;
    c = cw; /* 保存原始码字 */
    for (i = 1; i <= 12; i++)  /* 检查每个数据位 */
    {
        if (cw & 1)        /* 测试数据位 */
            cw ^= POLY;        /* 异或多项式 */
        cw >>= 1;            /* 移位中间结果 */
    }
    return ((cw << 12) | c);    /* 组装码字 */
}

int parity(unsigned long cw)
/* 此函数检查码字cw的整体奇偶校验。
   如果奇偶校验为偶数，则返回0，否则返回1。 */
{
    unsigned char p;

    /* 对码字的字节进行异或操作 */
    p = *(unsigned char *) &cw;
    p ^= *((unsigned char *) &cw + 1);
    p ^= *((unsigned char *) &cw + 2);

    /* 对中间结果的两半进行异或操作 */
    p = p ^ (p >> 4);
    p = p ^ (p >> 2);
    p = p ^ (p >> 1);

    /* 返回奇偶校验结果 */
    return (p & 1);
}

unsigned long syndrome(unsigned long cw)
/* 此函数计算并返回[23,12] Golay码字的综合。 */
{
    int i;
    cw &= 0x7fffffl;
    for (i = 1; i <= 12; i++)  /* 检查每个数据位 */
    {
        if (cw & 1)        /* 测试数据位 */
            cw ^= POLY;        /* 异或多项式 */
        cw >>= 1;            /* 移动中间结果 */
    }
    return (cw << 12);        /* 值与cw的高位配对 */
}

int weight(unsigned long cw)
/* 此函数计算23位码字cw的权重。 */
{
    int bits, k;

    /* Nibble权重表 */
    const char wgt[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

    bits = 0; /* 位计数器 */
    k = 0;
    /* 处理所有位，最多六个nibble */
    while ((k < 6) && (cw)) {
        bits = bits + wgt[cw & 0xf];
        cw >>= 4;
        k++;
    }

    return (bits);
}

unsigned long rotate_left(unsigned long cw, int n)
/* 此函数将23位码字cw向左旋转n位。 */
{
    int i;

    if (n != 0) {
        for (i = 1; i <= n; i++) {
            if ((cw & 0x400000l) != 0)
                cw = (cw << 1) | 1;
            else
                cw <<= 1;
        }
    }

    return (cw & 0x7fffffl);
}

unsigned long rotate_right(unsigned long cw, int n)
/* 此函数将23位码字cw向右旋转n位。 */
{
    int i;

    if (n != 0) {
        for (i = 1; i <= n; i++) {
            if ((cw & 1) != 0)
                cw = (cw >> 1) | 0x400000l;
            else
                cw >>= 1;
        }
    }

    return (cw & 0x7fffffl);
}

unsigned long correct(unsigned long cw, int *errs)
/* 此函数纠正Golay [23,12]码字cw，返回已纠正的码字。此函数将处理三个或更少的错误。对于四个或更多的错误，它将产生一些其他有效的Golay码字，可能不是预期的码字。*errs设置为纠正的位错误数。 */
{
    unsigned char
            w;                /* 当前综合限制权重，2或3 */
    unsigned long
            mask;             /* 用于翻转位的掩码 */
    int
            i, j;              /* 索引 */
    unsigned long
            s,                /* 计算的综合 */
    cwsaver;          /* 保存cw的初始值 */

    cwsaver = cw;         /* 保存 */
    *errs = 0;
    w = 3;                /* 初始综合权重阈值 */
    j = -1;               /* -1 = 第一次通过时没有尝试翻转位 */
    mask = 1;
    while (j < 23) /* 翻转每个试验位 */
    {
        if (j != -1) /* 切换一个试验位 */
        {
            if (j > 0) /* 恢复上一个试验位 */
            {
                cw = cwsaver ^ mask;
                mask += mask; /* 指向下一个位 */
            }
            cw = cwsaver ^ mask; /* 翻转下一个试验位 */
            w = 2; /* 在位操作时降低阈值 */
        }

        s = syndrome(cw); /* 查找错误 */
        if (s) /* 错误存在 */
        {
            for (i = 0; i < 23; i++) /* 检查每个循环移位的综合 */
            {
                if ((*errs = weight(s)) <= w) /* 综合匹配错误模式 */
                {
                    cw = cw ^ s;              /* 移除错误 */
                    cw = rotate_right(cw, i);  /* 反向旋转数据 */

                    if (j >= 0) /* 计算翻转的位数（按照Steve Duncan的说法） */
                        *errs = *errs + 1;

                    return (s = cw);
                } else {
                    cw = rotate_left(cw, 1);   /* 旋转到下一个模式 */
                    s = syndrome(cw);         /* 计算新的综合 */
                }
            }
            j++; /* 切换下一个试验位 */
        } else
            return (cw); /* 返回已纠正的码字 */
    }

    return (cwsaver); /* 如果没有纠正，则返回原始值 */
} /* correct */

int decode(int correct_mode, int *errs, unsigned long *cw)
/* 此函数以两种模式之一解码码字 *cw。如果correct_mode非零，将尝试纠正错误，*errs设置为纠正的位数，如果没有错误存在则返回0，如果存在奇偶校验错误则返回1。
如果correct_mode为零，则在 *cw 上执行错误检测，如果没有错误存在则返回0，如果存在总奇偶校验错误则返回1，如果存在码字错误则返回2。 */
{
    unsigned long parity_bit;

    if (correct_mode)               /* 纠正错误 */
    {
        parity_bit = *cw & 0x800000l; /* 保存奇偶校验位 */
        *cw &= ~0x800000l;            /* 移除奇偶校验位以进行纠正 */

        *cw = correct(*cw, errs);     /* 纠正最多三位错误 */
        *cw |= parity_bit;            /* 恢复奇偶校验位 */

        /* 检查是否有4位错误 */
        if (parity(*cw))            /* 奇偶校验错误 */
            return (1);
        return (0); /* 没有错误 */
    } else /* 仅检测错误 */
    {
        *errs = 0;
        if (parity(*cw)) /* 奇偶校验错误 */
        {
            *errs = 1;
            return (1);
        }
        if (syndrome(*cw)) {
            *errs = 1;
            return (2);
        } else
            return (0); /* 没有错误 */
    }
} /* 解码 */

void golay_test(void)
/* 这个函数测试了高雷(Golay)例程，用于检测和纠正各种错误限位的错误模式。error_mask循环遍历所有可能的值，并且error_limit选择了可以引起的最大错误数。*/
{
    unsigned long error_mask,        /* 用于引起错误的按位掩码 */
    trashed_codeword,  /* 试验校正的码字 */
    virgin_codeword;   /* 没有错误的原始码字 */
    unsigned char
            pass = 1,            /* 假设测试通过 */
    error_limit = 3;     /* 选择引起的位错误数量 */
    int error_count;       /* 接收纠正的错误数量 */

    virgin_codeword = golay(0x555); /* 制作一个测试码字 */
    if (parity(virgin_codeword))
        virgin_codeword ^= 0x800000l;
    for (error_mask = 0; error_mask < 0x800000l; error_mask++) {
        /* 过滤选定数量的位错误掩码 */
        if (weight(error_mask) <= error_limit) /* 你可以加快这个过程！ */
        {
            trashed_codeword = virgin_codeword ^ error_mask; /* 引起位错误 */

            decode(1, &error_count, &trashed_codeword); /* 尝试纠正位错误 */

            if (trashed_codeword ^ virgin_codeword) {
                printf("无法纠正由错误掩码=0x%lX引起的%d个错误\n", weight(error_mask), error_mask);
                pass = 0;
            }

//            if (kbhit()) /* 寻找用户输入 */
//            {
//                if (getch() == 27) return; /* 按下esc退出 */
//
//                /* 其他按键打印状态 */
//                printf("当前测试计数=%ld of %ld\n",error_mask,0x800000l);
//            }
        }
    }
    printf("Golay测试%s！\n", pass ? "通过" : "失败");
}


#endif //GOLAYCODE_GOGLAY_H
