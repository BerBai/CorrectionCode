/******************* main.cpp *******************/

#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <string>
#include "golay.h"
#include <ctime>

using namespace std;

int main() {
    system("chcp 65001");

    matrix gen,gen_t,par,par_t,encoded_message,syndrome,
            sB,s_add_colB,sB_add_colB_T,original_message;
    matrix message(1,12);
    matrix temp1(1,12);
    matrix noise(1,24);
    matrix error(1,24);
    matrix s_add_colB_full(12,13);
    matrix sB_add_colB_T_full(12,13);
    matrix temp(12,1);

    int i,j,k;
    vector<int> zeros(24);

    clock_t start = clock();

    gen = golay_generator();   // 生成 Golay 24 码的生成矩阵
    cout << "生成 Golay 24 码的生成矩阵（G）：" << endl;
    print(gen);

    gen_t = transpose(gen); // 生成矩阵 G 的转置
    cout << endl << "生成矩阵转置（G^T）：" << endl;
    print(gen_t);

    par = generate_B();        // Golay 24 码的奇偶校验矩阵
    cout << endl << "Golay 24 码的奇偶校验矩阵（B）：" << endl;
    print(par);

    par_t = transpose(par); // 奇偶校验矩阵 B 的转置
    cout << endl << "奇偶校验矩阵转置（B^T）：" << endl;
    print(par_t);

    message = generate_message();   // 生成随机消息
    cout << endl << "消息（m）：   ";
    print(message);

    encoded_message = encode(message,gen);  // 编码消息
    cout << endl << "编码消息（mG）：  ";
    print(encoded_message);

    noise = noisy_channel(encoded_message); // 将消息发送到“有噪声的信道”
    cout << endl << "接收到的消息（r = c + e）：  ";
    print(noise);

    syndrome = multiply(noise, gen_t);  // 计算 s = rG^T
    cout << endl << "综合（s = rG^T）：   ";
    print(syndrome);
    cout << "综合的权重是 " << weight(syndrome) << endl;

    sB = multiply(syndrome,par);    // 计算 sB
    cout << endl << "sB：    ";
    print(sB);
    cout  << "sB 的权重是 " << weight(sB);

    /* 计算行向量 s + (c_j)^T 并计算行向量的权重 */
    for(i = 0 ; i < 12 ; i++) {
        for(j = 0 ; j < 12 ; j++) {
            temp.m[j][0] = par.m[j][i];
        }

        temp1 = transpose(temp);
        s_add_colB = addition(syndrome,temp1);

        for(k = 0 ; k < 13 ; k++) {
            if(k < 12) {
                s_add_colB_full.m[i][k] = s_add_colB.m[0][k];
            } else {
                s_add_colB_full.m[i][k] = weight(s_add_colB);
            }
        }
    }

    cout << "\n\n" << "完整的 s + (c_j)^T（最后一列的值是行的权重）：" << endl;
    print(s_add_colB_full);

    /* 计算行向量 sB + (b_j)^T 并计算行向量的权重 */
    for(i = 0 ; i < 12 ; i++) {
        for(j = 0 ; j < 12 ; j++) {
            temp.m[j][0] = par_t.m[j][i];
        }

        temp1 = transpose(temp);
        sB_add_colB_T = addition(sB,temp1);

        for(k = 0 ; k < 13 ; k++) {
            if(k < 12) {
                sB_add_colB_T_full.m[i][k] = sB_add_colB_T.m[0][k];
            } else {
                sB_add_colB_T_full.m[i][k] = weight(sB_add_colB_T);
            }
        }
    }

    cout << "\n\n" << "完整的 sB + (b_j)^T（最后一列的值是行的权重）：" << endl;
    print(sB_add_colB_T_full);

    /***** 计算错误向量 e *****/
    if(weight(syndrome) <= 3) {
        for(i = 0; i < syndrome.m[0].size() ; i++) {
            if(syndrome.m[0][i] == 1) {
                error.m[0][i] = 1;
            }
        }
    } else if(weight(sB) <= 3) {
        for(i = 0; i < sB.m[0].size() ; i++) {
            if(sB.m[0][i] == 1) {
                error.m[0][i+12] = 1;
            }
        }
    }

    if(error.m[0] == zeros) {
        for(i = 0; i < 12 ; i++) {
            if(s_add_colB_full.m[i][12] <= 2) {
                error.m[0][i+12] = 1;

                for(j = 0 ; j< 12 ; j++) {
                    if(s_add_colB_full.m[i][j] ==1) {
                        error.m[0][j] = 1;
                    }
                }
            }
        }
    }

    if(error.m[0] == zeros) {
        for(i = 0; i < 12 ; i++) {
            if(sB_add_colB_T_full.m[i][12] <= 2) {
                error.m[0][i] = 1;

                for(j = 0 ; j< 12 ; j++) {
                    if(sB_add_colB_T_full.m[i][j] ==1) {
                        error.m[0][j+12] = 1;
                    }
                }
            }
        }
    }

    cout << endl <<"错误向量是：    ";
    print(error);

    original_message =addition(noise,error);    // 计算原始消息 c = r + e

    /***** 检查解码是否正确 *****/
    if(original_message.m == encoded_message.m) {
        cout << endl << "传输和纠错成功。" << endl;
    } else {
        cout << "传输或纠错出现错误。" << endl;
    }

    cout << endl << "计算出的原始发送消息：   ";
    print(original_message);
    cout << "         原始发送消息：   ";
    print(encoded_message);


    cout << endl << "计算出的原始消息 m ：   ";
    for(i = 0 ; i < 12 ; i++) {
        cout << original_message.m[0][i] << "  ";
    }
    cout << endl << "         原始消息 m ：   ";
    for(i = 0 ; i < 12 ; i++) {
        cout << message.m[0][i] << "  ";
    }

    clock_t finish = clock();

    cout << "\n\n" << "运行时间： " << (double)(finish-start)*1000/CLOCKS_PER_SEC<< " ms" << endl;

    return 0;
}
