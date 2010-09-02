/** @file rtkvarianobigeometry_ggo.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.3
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef RTKVARIANOBIGEOMETRY_GGO_H
#define RTKVARIANOBIGEOMETRY_GGO_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_PACKAGE "rtk"
#endif

#ifndef CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_PACKAGE_NAME "rtk"
#endif

#ifndef CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_RTKVARIANOBIGEOMETRY_VERSION "Creates an RTK geometry file from a Varian OBI acquisition."
#endif

/** @brief Where the command line options are stored */
struct args_info_rtkvarianobigeometry
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * config_arg;	/**< @brief Config file.  */
  char * config_orig;	/**< @brief Config file original value given at command line.  */
  const char *config_help; /**< @brief Config file help description.  */
  char * xml_file_arg;	/**< @brief Varian OBI XML information file on projections.  */
  char * xml_file_orig;	/**< @brief Varian OBI XML information file on projections original value given at command line.  */
  const char *xml_file_help; /**< @brief Varian OBI XML information file on projections help description.  */
  char * path_arg;	/**< @brief Path containing projections.  */
  char * path_orig;	/**< @brief Path containing projections original value given at command line.  */
  const char *path_help; /**< @brief Path containing projections help description.  */
  char * regexp_arg;	/**< @brief Regular expression to select projection files in path.  */
  char * regexp_orig;	/**< @brief Regular expression to select projection files in path original value given at command line.  */
  const char *regexp_help; /**< @brief Regular expression to select projection files in path help description.  */
  char * output_arg;	/**< @brief Output file name.  */
  char * output_orig;	/**< @brief Output file name original value given at command line.  */
  const char *output_help; /**< @brief Output file name help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int config_given ;	/**< @brief Whether config was given.  */
  unsigned int xml_file_given ;	/**< @brief Whether xml_file was given.  */
  unsigned int path_given ;	/**< @brief Whether path was given.  */
  unsigned int regexp_given ;	/**< @brief Whether regexp was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_rtkvarianobigeometry_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure args_info_rtkvarianobigeometry (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure args_info_rtkvarianobigeometry (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *args_info_rtkvarianobigeometry_purpose;
/** @brief the usage string of the program */
extern const char *args_info_rtkvarianobigeometry_usage;
/** @brief all the lines making the help output */
extern const char *args_info_rtkvarianobigeometry_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_rtkvarianobigeometry (int argc, char **argv,
  struct args_info_rtkvarianobigeometry *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_rtkvarianobigeometry_ext() instead
 */
int cmdline_parser_rtkvarianobigeometry2 (int argc, char **argv,
  struct args_info_rtkvarianobigeometry *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_rtkvarianobigeometry_ext (int argc, char **argv,
  struct args_info_rtkvarianobigeometry *args_info,
  struct cmdline_parser_rtkvarianobigeometry_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_rtkvarianobigeometry_dump(FILE *outfile,
  struct args_info_rtkvarianobigeometry *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_rtkvarianobigeometry_file_save(const char *filename,
  struct args_info_rtkvarianobigeometry *args_info);

/**
 * Print the help
 */
void cmdline_parser_rtkvarianobigeometry_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_rtkvarianobigeometry_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_rtkvarianobigeometry_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_rtkvarianobigeometry_params_init(struct cmdline_parser_rtkvarianobigeometry_params *params);

/**
 * Allocates dynamically a cmdline_parser_rtkvarianobigeometry_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_rtkvarianobigeometry_params structure
 */
struct cmdline_parser_rtkvarianobigeometry_params *cmdline_parser_rtkvarianobigeometry_params_create(void);

/**
 * Initializes the passed args_info_rtkvarianobigeometry structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_rtkvarianobigeometry_init (struct args_info_rtkvarianobigeometry *args_info);
/**
 * Deallocates the string fields of the args_info_rtkvarianobigeometry structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_rtkvarianobigeometry_free (struct args_info_rtkvarianobigeometry *args_info);

/**
 * The config file parser (deprecated version)
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_rtkvarianobigeometry_config_file() instead
 */
int cmdline_parser_rtkvarianobigeometry_configfile (const char *filename,
  struct args_info_rtkvarianobigeometry *args_info,
  int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_rtkvarianobigeometry_config_file (const char *filename,
  struct args_info_rtkvarianobigeometry *args_info,
  struct cmdline_parser_rtkvarianobigeometry_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_rtkvarianobigeometry_required (struct args_info_rtkvarianobigeometry *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* RTKVARIANOBIGEOMETRY_GGO_H */
