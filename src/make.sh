#!/bin/bash

ulimit -s  102400

cd ..
make
cd ./src



#./hydrocode.out Sod_test Sod_test/Sod_test_ROE 1 1_Roe free_1D
#./hydrocode.out Sod_10_test2 Sod_10_test2/Sod_10_test_ROE 2 1_Roe t1
#./hydrocode.out	Sod_10_test	Sod_10_test/Sod_10_test_ROE	2 1_Roe		Sod
#./hydrocode.out	Sod_10_test	Sod_10_test/Sod_10_test	2 2_GRP		Sod

#./hydrocode.out odd_even	odd_even/odd_even	2 1_Riemann_exact odd_even
#./hydrocode.out odd_even 	odd_even/odd_even_Roe 	2 1_Roe 		odd_even

#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 1_Riemann_exact	oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE_Roe	2 1_Roe			oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 2_GRP			oblique_periodic

#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_Roe	2 1_Roe			free
#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad	2 1_Riemann_exact	free
./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad	2 2_GRP			free
