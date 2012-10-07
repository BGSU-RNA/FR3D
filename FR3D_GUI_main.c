/*
 * MATLAB Compiler: 4.6 (R2007a)
 * Date: Mon Mar 19 14:19:34 2007
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe" "FR3D_GUI" 
 */

#include <stdio.h>
#include "mclmcr.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_FR3D_GUI_component_data;

#ifdef __cplusplus
}
#endif

static HMCRINSTANCE _mcr_inst = NULL;


#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_FR3D_GUI_C_API 
#define LIB_FR3D_GUI_C_API /* No special import/export declaration */
#endif

LIB_FR3D_GUI_C_API 
bool MW_CALL_CONV FR3D_GUIInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
    if (_mcr_inst != NULL)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!mclInitializeComponentInstance(&_mcr_inst,
                                        &__MCC_FR3D_GUI_component_data,
                                        true, NoObjectType, ExeTarget,
                                        error_handler, print_handler))
        return false;
    return true;
}

LIB_FR3D_GUI_C_API 
bool MW_CALL_CONV FR3D_GUIInitialize(void)
{
    return FR3D_GUIInitializeWithHandlers(mclDefaultErrorHandler,
                                          mclDefaultPrintHandler);
}

LIB_FR3D_GUI_C_API 
void MW_CALL_CONV FR3D_GUITerminate(void)
{
    if (_mcr_inst != NULL)
        mclTerminateInstance(&_mcr_inst);
}

int run_main(int argc, const char **argv)
{
    int _retval;
    /* Generate and populate the path_to_component. */
    char path_to_component[(PATH_MAX*2)+1];
    separatePathName(argv[0], path_to_component, (PATH_MAX*2)+1);
    __MCC_FR3D_GUI_component_data.path_to_component = path_to_component; 
    if (!FR3D_GUIInitialize()) {
        return -1;
    }
    _retval = mclMain(_mcr_inst, argc, argv, "FR3D_GUI", 1);
    if (_retval == 0 /* no error */) mclWaitForFiguresToDie(NULL);
    FR3D_GUITerminate();
    mclTerminateApplication();
    return _retval;
}

int main(int argc, const char **argv)
{
    if (!mclInitializeApplication(
        __MCC_FR3D_GUI_component_data.runtime_options,
        __MCC_FR3D_GUI_component_data.runtime_option_count))
        return 0;
    
    return mclRunMain(run_main, argc, argv);
}
