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

    initscr();
    cbreak();
    noecho();
    
    Screen screen(params.get_width(), params.get_height(), params.get_show_structure(), params.get_mode(), params.get_depthcharacter()); 
    screen.set_chainfile(params.get_chainfile(), params.get_in_file().size());
    for (int i = 0; i < params.get_in_file().size(); i++){
        screen.set_protein(params.get_in_file(i), i, params.get_show_structure());
    }
    screen.set_tmatrix();
    
    if (params.get_utmatrix() != ""){
        screen.set_utmatrix(params.get_utmatrix(),0);
    }
    screen.normalize_proteins(params.get_utmatrix());

    
    bool run = true;
    while(run) {
        screen.draw_screen();
        run = screen.handle_input();
        usleep(100);
    }

    endwin();
    return 0;
}