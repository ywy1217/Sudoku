close all;
% list of Sudokus (reproduced from Bradley Knockel)

% here is a typical easy Sudoku (one that can be solved by simply checking
% every row, column, box, and square over and over again) (I got this from
% the Yahoo game "Sudoku Daily")
matrix2=[
2,0,0,7,1,0,0,6,0
0,0,5,0,8,2,0,0,3
0,7,0,0,9,0,0,0,0
0,8,0,3,0,0,0,0,0
0,6,0,0,7,0,0,5,0
0,0,0,0,0,5,0,2,0
0,0,0,0,5,0,0,3,0
4,0,0,1,3,0,7,0,0
0,3,0,0,6,4,0,0,1];

% here is a moderate Sudoku (one that requires some fancy tricks) (some
% people have confused this one with the world's hardest)
matrix3=[
8,5,0,0,0,2,4,0,0
7,2,0,0,0,0,0,0,9
0,0,4,0,0,0,0,0,0
0,0,0,1,0,7,0,0,2
3,0,5,0,0,0,9,0,0
0,4,0,0,0,0,0,0,0
0,0,0,0,8,0,0,7,0
0,1,7,0,0,0,0,0,0
0,0,0,0,3,6,0,4,0];

% A hard Sudoku is one that requires tricks too fancy, too numerous,
% and too complex, so I resort to guessing. The following Sudoku is the
% famed "AI Escargot", which is claimed to be world's hardest, and I think
% it is.
matrix4=[
1,0,0,0,0,7,0,9,0
0,3,0,0,2,0,0,0,8
0,0,9,6,0,0,5,0,0
0,0,5,3,0,0,9,0,0
0,1,0,0,8,0,0,0,2
6,0,0,0,0,4,0,0,0
3,0,0,0,0,0,0,1,0
0,4,0,0,0,0,0,0,7
0,0,7,0,0,0,3,0,0];


puzzle = matrix4;

sudoku1 = Sudoku(puzzle);
% sudoku1.print_sudoku;
tic;
candidate_set = sudoku1.Sudoku_solver;
fprintf('my Algs (%.3fs)\n', toc);
candidate_set = candidate_set(:);
% sudoku1.print_sudoku(true);

tic;
ref_solution = sudoku1.sudoku_backtrack(puzzle);
fprintf('Ref Algs (%.3fs)\n', toc);

assert(all(sudoku1.solution == ref_solution,'all'));

