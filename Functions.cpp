//
// Created by vladu on 03.08.2021.
//

#include "Functions.h"

int save_state_file(array_2D temp_grid, const std::string& path){
    std::ofstream file(path);
    if(file){
        for(int i=0;i<temp_grid.size1();i++){
            for (int j = 0; j < temp_grid.size2(); ++j) {
                file<<temp_grid(i,j)<<" ";
            }
            file<<std::endl;
        }
        return 0;
    }
    else{
        return -1;
    }
}

Mpi_config::Mpi_config(const std::string& config, const std::string& starting_cond){
    std::fstream file(config);
    std::vector<std::string> text;
    if(file){
        std::string word;
        while( file >> word) {
            text.push_back(word);
        }
        file.close();
    }
    for(const std::string& line : text){
        std::string value_name = line.substr(0,line.find('='));
        if(value_name == "conductivity"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            conductivity = strtod(temp.c_str(), &end);
        }
        if(value_name == "height"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            height = strtod(temp.c_str(), &end);
        }
        if(value_name == "width"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            width = strtod(temp.c_str(), &end);
        }
        if(value_name == "delta_x"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            delta_x = strtod(temp.c_str(), &end);
        }
        if(value_name == "delta_y"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            delta_y = strtod(temp.c_str(), &end);
        }
        if(value_name == "delta_t"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            delta_t = strtod(temp.c_str(), &end);
        }
        if(value_name == "save_rate"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            save_rate = strtol(temp.c_str(), &end, 10);
        }
        if(value_name == "duration"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            cycle_duration = strtol(temp.c_str(), &end, 10);
        }
        if(value_name == "process_count"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            processes = strtol(temp.c_str(), &end,10);
        }
        if(value_name == "max_temperature"){
            auto temp = line.substr(line.find('='),line.size());
            char* end;
            temp.erase(0,1);
            max_temp = strtol(temp.c_str(), &end,10);
        }
    }
    std::fstream base_conditions(starting_cond);
    text.clear();
    if(base_conditions){
        std::string word;
        while( base_conditions >> word) {
            text.push_back(word);
        }
        base_conditions.close();
    }
    if(processes < 2){
        std::cerr<<"At least 2 processes are needed to run";
    }
    grid_height = round(height/delta_x);
    grid_width = round(width/delta_y);
    starting_condition.resize(grid_height, grid_width);
    int x = 0, y=0;

    for(const std::string& line: text){
        if(std::strtod(line.c_str(),nullptr) == 0 ||(x<1000 && y<1000)){
            starting_condition(x,y) = std::strtod(line.c_str(),nullptr);
            x++;
            if(x == grid_width-1){
                x=0;
                y++;
            }
        }
    }
}

/*
Magick::Color get_color_from_temp(double temperature, double max){
    Magick::Color output(255*(temperature/max),0,255*(1-temperature/max));
    return output;
}

int draw_image(twoD_array<double> temp_grid, const std::string& path, double max_temp){
    using namespace Magick;
    Image output(Magick::Geometry(temp_grid.h,temp_grid.w),Color(255, 255, 255, 0));
    for(int i=0;i<temp_grid.h;i++){
        for(int j=0;j<temp_grid.w;j++){
            output.pixelColor(i,j,get_color_from_temp(temp_grid.get(i,j), max_temp));
        }
    }
    output.write(path);
    return 0;
}
*/
