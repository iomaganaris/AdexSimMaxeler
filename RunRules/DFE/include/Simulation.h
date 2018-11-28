/**\file */
#ifndef SLIC_DECLARATIONS_Simulation_H
#define SLIC_DECLARATIONS_Simulation_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface readLMem ----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'readLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size (param_size * 4) bytes.
 */
void Simulation_readLMem(
	uint64_t param_size,
	uint64_t param_start,
	float *outstream_tocpu);

/**
 * \brief Basic static non-blocking function for the interface 'readLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size (param_size * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Simulation_readLMem_nonblock(
	uint64_t param_size,
	uint64_t param_start,
	float *outstream_tocpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'readLMem'
 * 
 */
typedef struct { 
	uint64_t param_size; /**<  [in] Interface Parameter "size". */
	uint64_t param_start; /**<  [in] Interface Parameter "start". */
	float *outstream_tocpu; /**<  [out] The stream should be of size (param_size * 4) bytes. */
} Simulation_readLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Simulation_readLMem_run(
	max_engine_t *engine,
	Simulation_readLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'readLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_readLMem_run_nonblock(
	max_engine_t *engine,
	Simulation_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Simulation_readLMem_run_group(max_group_t *group, Simulation_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_readLMem_run_group_nonblock(max_group_t *group, Simulation_readLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Simulation_readLMem_run_array(max_engarray_t *engarray, Simulation_readLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_readLMem_run_array_nonblock(max_engarray_t *engarray, Simulation_readLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Simulation_readLMem_convert(max_file_t *maxfile, Simulation_readLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*--------------------------- Interface writeLMem ----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'writeLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 4) bytes.
 */
void Simulation_writeLMem(
	uint64_t param_size,
	uint64_t param_start,
	const float *instream_fromcpu);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Simulation_writeLMem_nonblock(
	uint64_t param_size,
	uint64_t param_start,
	const float *instream_fromcpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMem'
 * 
 */
typedef struct { 
	uint64_t param_size; /**<  [in] Interface Parameter "size". */
	uint64_t param_start; /**<  [in] Interface Parameter "start". */
	const float *instream_fromcpu; /**<  [in] The stream should be of size (param_size * 4) bytes. */
} Simulation_writeLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Simulation_writeLMem_run(
	max_engine_t *engine,
	Simulation_writeLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMem_run_nonblock(
	max_engine_t *engine,
	Simulation_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Simulation_writeLMem_run_group(max_group_t *group, Simulation_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMem_run_group_nonblock(max_group_t *group, Simulation_writeLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Simulation_writeLMem_run_array(max_engarray_t *engarray, Simulation_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMem_run_array_nonblock(max_engarray_t *engarray, Simulation_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Simulation_writeLMem_convert(max_file_t *maxfile, Simulation_writeLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*-------------------------- Interface writeLMemInt --------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'writeLMemInt'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 1) bytes.
 */
void Simulation_writeLMemInt(
	uint64_t param_size,
	uint64_t param_start,
	const uint8_t *instream_fromcpu);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMemInt'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size (param_size * 1) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Simulation_writeLMemInt_nonblock(
	uint64_t param_size,
	uint64_t param_start,
	const uint8_t *instream_fromcpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMemInt'
 * 
 */
typedef struct { 
	uint64_t param_size; /**<  [in] Interface Parameter "size". */
	uint64_t param_start; /**<  [in] Interface Parameter "start". */
	const uint8_t *instream_fromcpu; /**<  [in] The stream should be of size (param_size * 1) bytes. */
} Simulation_writeLMemInt_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMemInt'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Simulation_writeLMemInt_run(
	max_engine_t *engine,
	Simulation_writeLMemInt_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMemInt'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMemInt_run_nonblock(
	max_engine_t *engine,
	Simulation_writeLMemInt_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMemInt'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Simulation_writeLMemInt_run_group(max_group_t *group, Simulation_writeLMemInt_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMemInt'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMemInt_run_group_nonblock(max_group_t *group, Simulation_writeLMemInt_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMemInt'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Simulation_writeLMemInt_run_array(max_engarray_t *engarray, Simulation_writeLMemInt_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMemInt'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_writeLMemInt_run_array_nonblock(max_engarray_t *engarray, Simulation_writeLMemInt_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Simulation_writeLMemInt_convert(max_file_t *maxfile, Simulation_writeLMemInt_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Auxiliary function to evaluate expression for "SimulationKernel.loopLength".
 */
int Simulation_get_SimulationKernel_loopLength( void );


/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_M_S Interface Parameter "M_S".
 * \param [in] param_M_T Interface Parameter "M_T".
 * \param [in] param_N Interface Parameter "N".
 * \param [in] param_N_adex Interface Parameter "N_adex".
 * \param [in] param_N_stdp Interface Parameter "N_stdp".
 * \param [in] param_Start Interface Parameter "Start".
 * \param [in] param_Steps Interface Parameter "Steps".
 * \param [in] param_WriteEn Interface Parameter "WriteEn".
 * \param [in] param_WriteStep Interface Parameter "WriteStep".
 * \param [in] param_burstSize Interface Parameter "burstSize".
 * \param [in] instream_adex_params The stream should be of size (param_N_adex * 4) bytes.
 * \param [in] instream_stdp_params The stream should be of size (param_N_stdp * 4) bytes.
 * \param [out] outstream_NeuronCPU The stream should be of size ((((param_WriteEn * param_M_T) * (param_Steps / param_WriteStep)) * 4) * 4) bytes.
 */
void Simulation(
	uint64_t param_M_S,
	uint64_t param_M_T,
	uint64_t param_N,
	int64_t param_N_adex,
	int64_t param_N_stdp,
	float param_Start,
	uint32_t param_Steps,
	uint32_t param_WriteEn,
	uint32_t param_WriteStep,
	int64_t param_burstSize,
	const float *instream_adex_params,
	const float *instream_stdp_params,
	float *outstream_NeuronCPU);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_M_S Interface Parameter "M_S".
 * \param [in] param_M_T Interface Parameter "M_T".
 * \param [in] param_N Interface Parameter "N".
 * \param [in] param_N_adex Interface Parameter "N_adex".
 * \param [in] param_N_stdp Interface Parameter "N_stdp".
 * \param [in] param_Start Interface Parameter "Start".
 * \param [in] param_Steps Interface Parameter "Steps".
 * \param [in] param_WriteEn Interface Parameter "WriteEn".
 * \param [in] param_WriteStep Interface Parameter "WriteStep".
 * \param [in] param_burstSize Interface Parameter "burstSize".
 * \param [in] instream_adex_params The stream should be of size (param_N_adex * 4) bytes.
 * \param [in] instream_stdp_params The stream should be of size (param_N_stdp * 4) bytes.
 * \param [out] outstream_NeuronCPU The stream should be of size ((((param_WriteEn * param_M_T) * (param_Steps / param_WriteStep)) * 4) * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Simulation_nonblock(
	uint64_t param_M_S,
	uint64_t param_M_T,
	uint64_t param_N,
	int64_t param_N_adex,
	int64_t param_N_stdp,
	float param_Start,
	uint32_t param_Steps,
	uint32_t param_WriteEn,
	uint32_t param_WriteStep,
	int64_t param_burstSize,
	const float *instream_adex_params,
	const float *instream_stdp_params,
	float *outstream_NeuronCPU);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t param_M_S; /**<  [in] Interface Parameter "M_S". */
	uint64_t param_M_T; /**<  [in] Interface Parameter "M_T". */
	uint64_t param_N; /**<  [in] Interface Parameter "N". */
	int64_t param_N_adex; /**<  [in] Interface Parameter "N_adex". */
	int64_t param_N_stdp; /**<  [in] Interface Parameter "N_stdp". */
	float param_Start; /**<  [in] Interface Parameter "Start". */
	uint32_t param_Steps; /**<  [in] Interface Parameter "Steps". */
	uint32_t param_WriteEn; /**<  [in] Interface Parameter "WriteEn". */
	uint32_t param_WriteStep; /**<  [in] Interface Parameter "WriteStep". */
	int64_t param_burstSize; /**<  [in] Interface Parameter "burstSize". */
	const float *instream_adex_params; /**<  [in] The stream should be of size (param_N_adex * 4) bytes. */
	const float *instream_stdp_params; /**<  [in] The stream should be of size (param_N_stdp * 4) bytes. */
	float *outstream_NeuronCPU; /**<  [out] The stream should be of size ((((param_WriteEn * param_M_T) * (param_Steps / param_WriteStep)) * 4) * 4) bytes. */
} Simulation_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Simulation_run(
	max_engine_t *engine,
	Simulation_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_run_nonblock(
	max_engine_t *engine,
	Simulation_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Simulation_run_group(max_group_t *group, Simulation_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_run_group_nonblock(max_group_t *group, Simulation_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Simulation_run_array(max_engarray_t *engarray, Simulation_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Simulation_run_array_nonblock(max_engarray_t *engarray, Simulation_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Simulation_convert(max_file_t *maxfile, Simulation_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Simulation_init(void);

/* Error handling functions */
int Simulation_has_errors(void);
const char* Simulation_get_errors(void);
void Simulation_clear_errors(void);
/* Free statically allocated maxfile data */
void Simulation_free(void);
/* These are dummy functions for hardware builds. */
int Simulation_simulator_start(void);
int Simulation_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Simulation_H */

