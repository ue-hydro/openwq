

#include "models_CH/headerfile_CH.hpp"

void OpenWQ_CH_model::CH_driver_config(
        OpenWQ_json& OpenWQ_json,
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_units& OpenWQ_units,
        OpenWQ_output& OpenWQ_output){
    
    // Local variable
    std::string msg_string;

    // NATIVE Bigoeochemical model
    // Parse biogeochemical expressions (and save in global)
    if ((OpenWQ_wqconfig.CH->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
    {
        
        bgc_flex_setBGCexpressions(
            OpenWQ_json,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_units,
            OpenWQ_output);

    }else{

        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: No BGC_module found or unkown";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
        
        // Abort (Fatal error)
        exit(EXIT_FAILURE);

    }
}