#pragma once

void Fill_x_trued();

void Fill_x_estd();

void Fill_u();

void Fill_C();

void Fill_Q();

void Fill_P();

void Fill_V();

Matrix* Calc_x_truec();

void GnuPlot_InitValues();

Matrix* Calc_x_estc();

Matrix* Calc_A_temp(Matrix* A);

Matrix* Calc_W();

void GnuPlotStuff1(int rep);

Matrix* Calc_A();

void print_alloc();