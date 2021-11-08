classdef Sudoku < handle
    properties
        puzzle
        solution
    end
    
    methods
        function obj = Sudoku(puzzle)
            obj.puzzle = puzzle;
        end
    end
    
    methods(Static)
        solution = sudoku_backtrack(puzzle)
    end
end