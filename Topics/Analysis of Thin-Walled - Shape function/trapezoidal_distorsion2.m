clear;
clc;
syms a b h sina cosa real

k=[(a+b)/4 a/2 (a+b)/4 b/2
    h*cosa/sina a -h*cosa/sina -b
    1 0 1 0
  1 0 0 0];
  f=[0 0 0 1]';
d1=k\f
