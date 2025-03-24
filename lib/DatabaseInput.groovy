// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.

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
        throw new RuntimeException(message)
    }

    static void assertDir(String dir) {
        def path = new File(dir)
        if (!path.exists() || !path.isDirectory()) {
            exitError("The directory '${path}' does not exist. Exiting.")
        }
    }

    static void assertFile(String fname) {
        def path = new File(fname)
        if (!path.exists() || !path.isFile()) {
            exitError("The file '${path}' does not exist. Exiting.")
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
            exitError("Missing required key '$key' in configuration")
        }
        def fname = "${path}/${config[key]}"
        assertFile(fname)
        return fname
    }

    static void checkUniqueKeys(List<Map<String, Object>> maps) {
        def seenKeys = new HashSet<String>()
        maps.each { map ->
            map.each { key, value ->
                def normalizedKey = key.toLowerCase()
                if (seenKeys.contains(normalizedKey)) {
                    exitError("Duplicate key detected: $key")
                }
                seenKeys.add(normalizedKey)
            }
        }
    }

    static Map upperCaseMap(List<Map<String, Object>> maps) {
        checkUniqueKeys(maps)
        def merged = maps.inject([:]) { mergedSoFar, map -> mergedSoFar + map }
        def mapWithUpperCaseKeys = merged.collectEntries {
            name, fasta -> [(name.toUpperCase()): fasta]
        }
        return mapWithUpperCaseKeys
    }

    static Map getVirusFilterPaths(Map yamlConfig) {
        def all = ["ALL": yamlConfig.all.fasta]
        def family_filters = yamlConfig.filters
        def curated = yamlConfig.curated.collectEntries {
            name, entries -> [(name): entries.fasta]
        }
        return upperCaseMap([all, family_filters, curated])
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
        def contaminationConfig = upperCaseMap([readYamlConfig(contaminationConfigFileName)])
        this.contaminationFilters = dbParams.contamination_filters.tokenize(",")
        this.contaminationFilterFiles = this.contaminationFilters.collect {
            contaminant -> getFileFromConfig(contaminationConfig, this.contaminantsDir, contaminant.toUpperCase())
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
            target -> [
                target: target,
                path: getFileFromConfig(filterFilenames, virusDir, target.toUpperCase())
            ]
        }

        // blast
        this.blastDir = "${this.virusDir}/${virusConfig.all.blast_db}"
        assertDir(this.blastDir)

        this.blastPrefix = virusConfig.all.blast_prefix
        assertFile("${this.blastDir}/${this.blastPrefix}.ndb")

        // classification
        this.classificationLibrary = dbParams.centrifuge_classification_library
        assertFile("${this.classificationDir}/${this.classificationLibrary}.1.cf")
    }
}
