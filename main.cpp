#include "Functions.h"
#include <boost/mpi.hpp>

void get_col(double* out, array_2D input, int num){
    for(int j=0; j<input.size2();j++){
        out[j] = input(num,j);
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Mpi_config conf;
    if(argc > 3){
        conf = Mpi_config(argv[1],argv[2]);
    }
    else if(argc == 2){
        conf = Mpi_config(argv[1]);
    }
    else{
        conf = Mpi_config();
    }

    if (rank == 0) {
        int width = conf.grid_width/(conf.processes -1);
        array_2D grid(conf.grid_height,width);
        long long cycle = 0;
        double out[conf.grid_height];
        for(int i =0;i<conf.grid_height;i++){
            for(int j=0;j<width;j++){
                grid(i,j) = conf.starting_condition(i,j);
            }
        }
        std::cout<<conf.starting_condition.size1()<<" "<<conf.starting_condition.size2()<<std::endl;
        std::cout<<conf.starting_condition(0,999)<<std::endl;
        while(cycle < conf.cycle_duration){
            if(cycle%conf.save_rate == 0){
                for(int j=0; j<grid.size2();j++){
                    get_col(out, grid,j);
                    MPI_Send(out, conf.grid_height, MPI_DOUBLE, conf.processes-1, 0, MPI_COMM_WORLD);
                }
                bool done = true;
                MPI_Send(&done, 1, MPI_C_BOOL, rank+1,0,MPI_COMM_WORLD);
            }
            for(int i =1; i<conf.grid_height-1;i++){
                for(int j=1; j<width-1;j++){
                    double value =grid(i,j);
                    value += (grid(i+1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                    value += (grid(i-1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                    value += (grid(i,j+1)- value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                    value += (grid(i,j-1)- value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                    grid(i,j) = value;
                }
                int j = width - 1;
                double value =grid(i,j);
                double recv , send;
                MPI_Recv(&recv,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                value += (grid(i+1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                value += (grid(i-1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                value += (grid(i,j-1)- value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                value += (recv - value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                grid(i,j) = value;
                send= grid(i,j);
                MPI_Send(&send,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD);
            }
            cycle++;
        }
    } else if (rank == conf.processes-2){
        int width = conf.grid_width / (conf.processes-1);
        array_2D grid(conf.grid_height,width);
        long long cycle = 0;
        double out[conf.grid_height];
        bool done;
        for (int i = 0; i < conf.grid_height; i++) {
            for (int j = 0; j < width; j++) {
                grid(i,j) = conf.starting_condition(i,j+(conf.processes-2)*(conf.grid_width/(conf.processes-1)));
            }
        }
        while (cycle < conf.cycle_duration) {
            if (cycle%conf.save_rate == 0) {
                MPI_Recv(&done, 1, MPI_C_BOOL, rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int j=0; j<grid.size2();j++){
                    get_col(out, grid, j);
                    MPI_Send(out, conf.grid_height, MPI_DOUBLE, conf.processes-1, 0, MPI_COMM_WORLD);
                }
            }
            for (int i = 1; i < conf.grid_height - 1; i++) {
                int j = 0;
                double value =grid(i,j);
                double recv, send;
                send= grid(i,j);
                MPI_Send(&send,1,MPI_DOUBLE,conf.processes-3,1,MPI_COMM_WORLD);
                MPI_Recv(&recv,1,MPI_DOUBLE,conf.processes-3,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                value += (grid(i+1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                value += (grid(i-1,j)- value)*conf.conductivity*conf.delta_t*conf.delta_x*conf.delta_x;
                value += (grid(i,j+1)- value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                value += (recv - value)*conf.conductivity*conf.delta_t*conf.delta_y*conf.delta_y;
                grid(i,j) = value;
                for (j = 1; j < width-1; j++) {
                    value = grid(i,j);
                    value += (grid(i + 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                    value += (grid(i - 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                    value += (grid(i,j + 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y*conf.delta_y;
                    value += (grid(i,j - 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y*conf.delta_y;
                    grid(i,j) = value;
                }
            }
            cycle++;
        }
    } else if (rank == conf.processes-1){
        //saver process
        double in[conf.grid_height];
        long long iterations = conf.cycle_duration/conf.save_rate, cur_iter =0;
        while(cur_iter< iterations) {
            std::cout<<"Current iteration: "<<cur_iter<<std::endl;
            std::vector<array_2D> matrices;
            bool transferred = true;
            array_2D matrix(conf.grid_height,conf.grid_width);
            matrices.reserve(conf.processes-1);
            for(int i=0;i<conf.processes-1;i++){
                matrices.emplace_back(conf.grid_height,conf.grid_width/(conf.processes-1));
            }
            for(int i=0;i<conf.processes-1;i++){
                    for(int k=0;k<matrices[i].size2();k++){
                        if(!MPI_Recv(in, conf.grid_height, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE)){
                            transferred = false;
                        }
                        for(int j=0; j< conf.grid_width/(conf.processes-1);j++)
                        matrices[i](j,k) = in[j];
                    }
                std::cout<<"Data from process "<<i<<" saved"<<std::endl;
            }
            std::cout<<"Transfer status: "<<transferred<<std::endl;
                for (int k = 0; k < conf.processes-2; k++) {
                    for (int i = 0; i < conf.grid_height; i++) {
                        for (int j = 0; j < conf.grid_width / (conf.processes-1); j++) {
                            matrix(i,j + k * (conf.grid_width / (conf.processes-1))) =  matrices[k](i,j);
                        }
                    }
                //draw_image(matrix,"results/result"+std::to_string(cur_iter)+".png", max_temp);
                save_state_file(matrix,"results/result"+std::to_string(cur_iter)+".txt");
            }
            cur_iter += 1;
        }
    }
    else {
        array_2D grid(conf.grid_height,conf.grid_width / (conf.processes - 1));
        long long cycle = 0;
        double out[conf.grid_height];
        bool done;
        int width = conf.grid_width / (conf.processes - 1);
        for (int i = 0; i < conf.grid_height; i++) {
            for (int j = 0; j < width; j++) {
                grid(i,j) = conf.starting_condition(i,j + rank * (conf.grid_width / (conf.processes - 1)));
            }
        }
        while (cycle < conf.cycle_duration) {
            if (cycle%conf.save_rate == 0) {
                MPI_Recv(&done, 1, MPI_C_BOOL, rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int j=0; j<grid.size2();j++){
                    get_col(out, grid, j);
                    MPI_Send(out, conf.grid_height, MPI_DOUBLE, conf.processes-1, 0, MPI_COMM_WORLD);
                }
                MPI_Send(&done, 1, MPI_C_BOOL, rank+1,0,MPI_COMM_WORLD);
            }
            for (int i = 1; i < conf.grid_height - 1; i++) {
                int j = 0;
                double value = grid(i,j);
                double recv, send;
                send = grid(i,j);
                MPI_Send(&send, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
                MPI_Recv(&recv, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                value += (grid(i + 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                value += (grid(i - 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                value += (grid(i,j + 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                value += (recv - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                grid(i,j) = value;
                for (j = 1; j < width - 1; j++) {
                    value = grid(i,j);
                    value += (grid(i + 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                    value += (grid(i - 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                    value += (grid(i,j + 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                    value += (grid(i,j - 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                    grid(i,j) = value;
                }
                j = width - 1;
                value = grid(i,j);
                MPI_Recv(&recv, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                value += (grid(i + 1,j) - value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                value += (grid(i - 1,j) -value) * conf.conductivity * conf.delta_t * conf.delta_x* conf.delta_x;
                value += (grid(i,j - 1) - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                value += (recv - value) * conf.conductivity * conf.delta_t * conf.delta_y* conf.delta_y;
                grid(i,j) = value;
                send = grid(i,j);
                MPI_Send(&send, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            }
            cycle++;
        }
    }
    std::cout<<"Process "<<rank<<" finished"<<std::endl;
    MPI_Finalize();
    return 0;
}


