package org.molgenis.asterix.pipeline;

import org.apache.commons.cli.ParseException;
import org.molgenis.asterix.config.ConfigProvider;

import java.io.IOException;

/**
 * Class that is the application entrypoint
 *
 * @author OelenR
 */
public class AppStarter {

    /**
     * Application entry
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        AppStarter app = new AppStarter();
        app.start(args);
    }

    /**
     * Start the application
     * @param args the command line arguments
     */
    public void start(String[] args) {
        try {
            //load config
            ConfigProvider configProvider = ConfigProvider.getInstance();
            configProvider.loadCliArguments(args);
            //check for help request by user
            if(!configProvider.isRequestedHelp()) {
                Pipeline pipeline = new Pipeline();
                pipeline.doPipeline();
            }
        }
        catch (IOException | ParseException e) {
            e.printStackTrace();
        }
    }

}
