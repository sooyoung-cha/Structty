#include <iostream>
#include <ncurses.h>
#include <unistd.h> 
#include "Protein.hpp"
#include "Parameters.hpp"
#include "Screen.hpp"

int main(int argc, char* argv[]) {
    Parameters params(argc, argv);

    if (!params.check_arg_okay()) {
        return -1; 
    }
    params.print_args();

    
    Screen screen(params.get_width(), params.get_height(), params.get_show_structure(), params.get_mode(), params.get_depthcharacter()); 
    
    for (int i = 0; i < params.get_in_file().size(); i++){
        screen.set_protein(params.get_in_file(i), params.get_chains(i), params.get_show_structure());
    }
    screen.set_tmatrix();
    
    if (params.get_utmatrix() != ""){
        screen.set_utmatrix(params.get_utmatrix(),0);
    }
    screen.normalize_proteins(params.get_utmatrix());

    initscr(); // ncurses 초기화
    cbreak();  // 특수 키 입력을 받도록 설정
    noecho();  // 입력한 문자를 화면에 표시하지 않음
    
    bool run = true;
    while(run) {
        screen.draw_screen();
        run = screen.handle_input();
        usleep(100);
    }

    endwin(); // ncurses 종료
    return 0;
}