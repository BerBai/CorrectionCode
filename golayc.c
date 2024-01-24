/*
golay.c: 实现二进制Golay码[24,12,8]。
改编自：http://aqdi.com/articles/using-the-golay-error-detection-and-correction-code-3/

作者：Axel Persinger
许可证：MIT许可证
*/

#include "golayc.h"

int add(int i, int j)
{
    return i + j;
}

/**
 * @brief  计算并填写码字的校验位
 * @note   这里是填写码字的校验位
 * @param  code: 码字信息
 * @retval None
 */
void calculate_check(st_codeword* code)
{
    code->check = 0;  // 算法开始时，检验位为0
    short c = code->information;
}

/**
 * @brief  计算并填写码字的奇偶校验位
 * @note   这里是填写码字的奇偶校验位
 * @param  code: 码字信息
 * @retval None
 */
void calculate_parity(st_codeword* code)
{
    code->parity = 0;  // 算法开始时，奇偶校验位为0
    short c = code->information;
    do
    {
        if (c & 0x1)
            code->parity = ~code->parity;  // 如果信息中有1，则奇偶校验位取反

    } while (c >>= 1);
}

/**
 * @brief  计算并返回码字的权重
 * @note   这里是计算并返回码字的权重
 * @param  code: 码字信息
 * @retval weight
 */
short calculate_weight(st_codeword* code)
{
    short weight = 0;
    short c = code->information;
    do
    {
        if (c & 0x1)
            weight += 1;  // 如果信息中有1，则权重加1

    } while (c >>= 1);

    return weight;
}

int main(int argc, char const *argv[])
{
    st_codeword* code;
    code = (st_codeword*) malloc(sizeof(st_codeword));

    code->information = 0b01101;  // 设置码字信息

    // printf("code->parity = %d\n", code->check);
    calculate_parity(code);  // 计算奇偶校验位
    // int i = 0x3;
    printf("weight(%d) = %d\n", code->information, calculate_weight(code));  // 打印权重

    return 0;
}
