include { kraken_db } from '../modules/kraken'

workflow download_kraken_db {
    main:
        if (params.kraken) {
            kraken_db_preload = file("${params.databases}/kraken/GRCh38.p13_SC2_2022-03-01")
            if (kraken_db_preload.exists()) { database_kraken = kraken_db_preload }    
            else { kraken_db() ; database_kraken = kraken_db.out }
        }
    emit: database_kraken
}
