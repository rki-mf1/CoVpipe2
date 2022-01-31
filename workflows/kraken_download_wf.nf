include { kraken_db } from '../modules/kraken'

workflow download_kraken_db {
    main:
        if (params.kraken) {
            // local storage via storeDir
            if (!params.cloudProcess) { kraken_db() ; database_kraken = kraken_db.out}
            // cloud storage file.exists()?
            if (params.cloudProcess) { 
                kraken_db_preload = file("${params.databases}/GRCh38.p13_SC2_2021-02-08")
                if (kraken_db_preload.exists()) { database_kraken = kraken_db_preload }    
                else { kraken_db() ; database_kraken = kraken_db.out }
            }
        }
    emit: database_kraken
}
