import org.yaml.snakeyaml.Yaml


class DatabaseInput {

    String virusDir
    String virusConfigFileName
    String contaminantsDir
    String classificationDir
    String blastDir
    String blastPrefix
    List<String> contaminationFilters
    List<String> contaminationFilterFiles
    List<Map<String, String>> virusTargets
    List<String> classificationLibraries

    static void exitError(String message) {
        System.err.println(message)
        System.exit(1)
    }

    static void assertDir(String dir) {
        def path = new File(dir)
        if (!path.exists() || !path.isDirectory()) {
            exitError("ERROR: The directory '${path}' does not exist. Exiting.")
        }
    }

    static void assertFile(String fname) {
        def path = new File(fname)
        if (!path.exists() || !path.isFile()) {
            exitError("ERROR: The file '${path}' does not exist. Exiting.")
        }
    }

    static String getDir(String dir, String defaultDir) {
        def outdir = dir ?: defaultDir
        assertDir(outdir)
        return outdir
    }

    static String getFile(String fname, String defaultFname) {
        def fnameOut = fname ?: defaultFname
        assertFile(fnameOut)
        return fnameOut
    }

    static Map readYamlConfig(String fname) {
        def yamlFile = new File(fname)
        def yamlParser = new Yaml()
        def configData = yamlParser.load(yamlFile.text)
        return configData
    }

    static String getFileFromConfig(Map config, String path, String key) {
        if (!config.containsKey(key)) {
            exitError("Error: Missing required key '$key' in configuration")
        }
        def fname = "${path}/${config[key]}"
        assertFile(fname)
        return fname
    }

    static Map getVirusFilterPaths(Map yamlConfig) {
        def all = [ALL: yamlConfig.all.fasta]
        def family_filters = yamlConfig.filters
        def curated = yamlConfig.curated.collectEntries { name, entries ->
            [(name): entries.fasta]
        }
        return all + family_filters + curated
    }

    DatabaseInput(Map dbParams) {
        def baseDir = getDir(dbParams.base_db, dbParams.database_defaults.base)
        
        this.virusDir = getDir(
            dbParams.virus_db,
            "${baseDir}/${dbParams.database_defaults.virus}"
        )
        this.contaminantsDir = getDir(
            dbParams.contaminants_db,
            "${baseDir}/${dbParams.database_defaults.contaminants}"
        )
        this.classificationDir = getDir(
            dbParams.classification_db,
            "${baseDir}/${dbParams.database_defaults.classification}"
        )

        // contamination
        def contaminationConfigFileName = getFile(
            dbParams.contaminants_db_config,
            "${this.contaminantsDir}/${dbParams.database_defaults.contaminants_db_config}"
        )
        def contaminationConfig = readYamlConfig(contaminationConfigFileName)
        this.contaminationFilters = dbParams.contamination_filters.tokenize(",")
        this.contaminationFilterFiles = this.contaminationFilters.collect {
            contaminant -> getFileFromConfig(contaminationConfig, this.contaminantsDir, contaminant)
        }

        // virus
        this.virusConfigFileName = getFile(
            dbParams.virus_db_config,
            "${this.virusDir}/${dbParams.database_defaults.virus_db_config}"
        )
        def virusConfig = readYamlConfig(this.virusConfigFileName)
        def filterFilenames = getVirusFilterPaths(virusConfig)
        def virusTargetNames = dbParams.targets.tokenize(",")
        this.virusTargets = virusTargetNames.collect {
            target -> [target: target, path: getFileFromConfig(filterFilenames, virusDir, target)]
        }

        // blast
        this.blastDir = "${this.virusDir}/${virusConfig.all.blast_db}"
        assertDir(this.blastDir)

        this.blastPrefix = virusConfig.all.blast_prefix
        assertFile("${this.blastDir}/${this.blastPrefix}.ndb")

        // classification
        this.classificationLibraries = dbParams.centrifuge_classification_libraries.tokenize(",")
    }
}
