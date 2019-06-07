package org.molgenis.asterix.org.molgenis.asterix.config;

import java.util.*;

import org.apache.commons.cli.*;

/**
 * @author OelenR
 *
 * class to store the configuration
 */
public class ConfigProvider {

    /*
    this class implements the singleton pattern
     */

    //loaded and default properties for configuration
    private Properties storedConfig = null;
    //parser for CLI arguments
    private CommandLineParser parser = null;
    //values to check for CLI args
    private List<String> possibleCliArguments = null;
    //mapping of cli argument to option name
    private Map<String, String> cliOptionToName;
    //options that are possible
    private Options options = null;


    //self reference
    private static ConfigProvider SELF = new ConfigProvider();

    //private constructor according to singleton pattern
    private ConfigProvider(){
        this.storedConfig = new Properties();
        this.parser = new DefaultParser();
        this.possibleCliArguments = new ArrayList<>();
        this.cliOptionToName = new HashMap<>();
        this.options = new Options();
        this.setPossibleOptions();
        this.setCreateOptions();
        this.setMapCliOptionsToName();
        this.setDefaultProperties();
    }

    /**
     * get an instance of the ConfigProvider
     * @return the ConfigProvider
     */
    public static ConfigProvider getInstance(){
        return SELF;
    }

    /**
     * add possible options to options list
     */
    private void setPossibleOptions(){
        this.possibleCliArguments.add(ConfigConstants.OPTION_STAR_ALLELE_OUTPUT_DIR);
        this.possibleCliArguments.add(ConfigConstants.OPTION_SNP_HAPLO_TABLE_DIR);
        this.possibleCliArguments.add(ConfigConstants.OPTION_HAPLOTYPE_DIR);
        this.possibleCliArguments.add(ConfigConstants.OPTION_HAPLO_PHENO_TABLE_DIR);
        this.possibleCliArguments.add(ConfigConstants.OPTION_PREDICTED_PHENOTYPES_OUTPUT_DIR);
    }

    /**
     * set options with descriptions
     */
    private void setCreateOptions(){
        this.options.addOption(ConfigConstants.OPTION_STAR_ALLELE_OUTPUT_DIR,true,"the output directory for the star alleles files");
        this.options.addOption(ConfigConstants.OPTION_SNP_HAPLO_TABLE_DIR,true,"the directory of the tables to map SNPs to haplotyes");
        this.options.addOption(ConfigConstants.OPTION_HAPLO_PHENO_TABLE_DIR, true, "the directory of the tables to map the haplotypes to phenotypes");
        this.options.addOption(ConfigConstants.OPTION_HAPLOTYPE_DIR, true, "the directory of the haplotypes to parse");
        this.options.addOption(ConfigConstants.OPTION_PREDICTED_PHENOTYPES_OUTPUT_DIR, true, "the directory where the predicted phenotypes are written");
    }

    /**
     * map CLI option strings to constant names
     */
    private void setMapCliOptionsToName(){
        this.cliOptionToName.put(ConfigConstants.OPTION_STAR_ALLELE_OUTPUT_DIR, ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        this.cliOptionToName.put(ConfigConstants.OPTION_SNP_HAPLO_TABLE_DIR, ConfigConstants.SNP_HAPLO_TABLE_DIR);
        this.cliOptionToName.put(ConfigConstants.OPTION_HAPLO_PHENO_TABLE_DIR, ConfigConstants.HAPLO_PHENO_TABLE_DIR);
        this.cliOptionToName.put(ConfigConstants.OPTION_HAPLOTYPE_DIR, ConfigConstants.HAPLOTYPE_DIR);
        this.cliOptionToName.put(ConfigConstants.OPTION_PREDICTED_PHENOTYPES_OUTPUT_DIR, ConfigConstants.PREDICTED_PHENOTYPES_OUTPUT_DIR);
    }

    /**
     * set the default properties
     */
    private void setDefaultProperties(){
        this.storedConfig.put(ConfigConstants.STAR_ALLELE_OUTPUT_DIR,ConfigConstants.DEFAULT_STAR_ALLELE_OUTPUT_DIR);
        this.storedConfig.put(ConfigConstants.SNP_HAPLO_TABLE_DIR, ConfigConstants.DEFAULT_SNP_HAPLO_TABLE_DIR);
        this.storedConfig.put(ConfigConstants.HAPLOTYPE_DIR, ConfigConstants.DEFAULT_HAPLOTYPE_DIR);
        this.storedConfig.put(ConfigConstants.HAPLO_PHENO_TABLE_DIR, ConfigConstants.DEFAULT_HAPLO_PHENO_TABLE_DIR);
        this.storedConfig.put(ConfigConstants.PREDICTED_PHENOTYPES_OUTPUT_DIR, ConfigConstants.DEFAULT_PREDICTED_PHENOTYPES_OUTPUT_DIR);
    }

    /**
     * load the command line arguments, overwriting the defaults
     * @param args the command line arguments
     */
    public void loadCliArguments(String[] args) throws ParseException {
        //null safe operation
        if(null != args && args.length>0){
            //parse options
            CommandLine cmd = parser.parse( options, args);
            //if properties are available, overwrite the defaults
            for(String option : this.possibleCliArguments){
                if(cmd.hasOption(option)){
                    //get constant name
                    String name = this.cliOptionToName.get(option);
                    //get value
                    String value = cmd.getOptionValue(option);
                    this.storedConfig.put(name, value);
                }
            }
        }
    }

    /**
     * check if parameter is loaded
     * @param paramName the name of the parameter
     * @return whether the parameter is loaded
     */
    public boolean containsConfigParam(String paramName){
        if (this.storedConfig.contains(paramName)) {
            return true;
        }
        else{
            return false;
        }
    }

    /**
     * get parameter value
     * @param paramName the name of the parameter
     * @return the value of the parameter, or null if not present
     */
    public String getConfigParam(String paramName){
        return this.storedConfig.getProperty(paramName);
    }

}
