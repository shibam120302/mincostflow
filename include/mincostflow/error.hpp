#pragma once

#include <stdexcept>
#include <string>

namespace ln
{
    enum error_codes { NO_ERROR=0, UNKNOWN_ERROR, UNFEASIBLE_FLOW, ARC_DUPLICATE, ARC_MISSING };
    
    // general purpose exception
    class mcf_exception : public std::exception
    {
        const std::string what_str{};
        
        public:
        mcf_exception()
        {}
        
        mcf_exception(const std::string& arg):
            what_str{arg}
        {}
        
        virtual int error_code()const = 0;
        
        const char* what()const noexcept override
        {
            return what_str.c_str();
        }
    };
    
    // exception to signal that no feasible flow could be found
    class mcf_unfeasible : public mcf_exception
    {
        static constexpr int my_error_code = error_codes::UNFEASIBLE_FLOW;
        public:
        mcf_unfeasible(const std::string& what_arg):
            mcf_exception{"Unfeasible flow, " + what_arg}
        {
        }
        mcf_unfeasible(){}
        
        int error_code()const override
        {
            return my_error_code;
        }
    };
    
    // exception to signal arc duplication
    class mcf_arc_duplicate : public mcf_exception
    {
        static constexpr int my_error_code = error_codes::ARC_DUPLICATE;
        public:
        mcf_arc_duplicate(const std::string& what_arg):
            mcf_exception{"Arc duplication, "+what_arg}
        {
        }
        mcf_arc_duplicate(){}
        
        int error_code()const override
        {
            return my_error_code;
        }
    };
    
    // exception to signal a request for a missing arc
    class mcf_arc_missing : public mcf_exception
    {
        static constexpr int my_error_code = error_codes::ARC_MISSING;
        public:
        mcf_arc_missing(const std::string& what_arg):
            mcf_exception{"Arc missing, "+what_arg}
        {
        }
        mcf_arc_missing(){}
        
        int error_code()const override
        {
            return my_error_code;
        }
    };
}
