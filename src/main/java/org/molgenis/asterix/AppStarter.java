package org.molgenis.asterix;

import org.apache.commons.cli.ParseException;
import org.molgenis.asterix.org.molgenis.asterix.config.ConfigProvider;

import java.io.IOException;

/**
 * class that is the application entrypoint
 */
public class AppStarter {

    /**
     * application entry
     * @param args the command line arguments
     */
    public static void main(String[] args){
        AppStarter app = new AppStarter();
        app.start(args);
    }

    /**
     * start the application
     * @param args the command line arguments
     */
    public void start(String[] args){
        //load config
        this.loadCliArgs(args);
        //do the actual work
        try {
            StarAlleleToPhenotype.run();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * load the command line arguments
     * @param args the command line arguments
     */
    private void loadCliArgs(String[] args){
        try {
            ConfigProvider configProvider = ConfigProvider.getInstance();
            configProvider.loadCliArguments(args);
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
    }
}
