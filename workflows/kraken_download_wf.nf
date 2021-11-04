include { kraken_db } from '../modules/kraken'

workflow download_kraken_db {
    main:
        if (params.kraken) {
            // local storage via storeDir
            if (!params.cloudProcess) { kraken_db() ; database_kraken = kraken_db.out}
            // cloud storage file.exists()?
            if (params.cloudProcess) { 
                kraken_db_preload = file("${params.databases}/kraken/GRCh38.p13_GBcovid19-2020-05-22.tar.gz")
                if (kraken_db_preload.exists()) { database_kraken = kraken_db_preload }    
                else { kraken_db() ; database_kraken = kraken_db.out }
            }
        }
    emit: database_kraken
}
