add_library(fabm_models_mops OBJECT
            mops_model_library.F90
            shared.F90
            carbon.F90
            oxygen.F90
            insolation.F90
            radiation.F90
            phytoplankton.F90
            zooplankton.F90
            remineralization.F90
            detritus.F90
            tracer.F90
            nitrogen_fixation.F90
            runoff.F90
           )

add_dependencies(fabm_models_mops fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2023 - BB
#-----------------------------------------------------------------------
