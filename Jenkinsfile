pipeline {
    agent { label 'monarch-agent-medium' }
        environment {
        HOME = "${env.WORKSPACE}"
        RELEASE = sh(script: "echo `date +%Y-%m-%d`", returnStdout: true).trim()
        PATH = "/opt/poetry/bin:${env.PATH}"
    }
    stages {
        stage('setup') {
            steps {
                sh 'poetry install'
            }
        }
        stage('generate-mapping-file') {
            steps {
                sh 'poetry run gene-mapping generate --download --preprocess-uniprot'
            }
        }
        stage('upload-mapping-file'){
            steps{
                sh '''
                    gsutil cp output/gene_mappings.tsv gs://monarch-archive/monarch-gene-mapping/${RELEASE}/
                    gsutil cp output/gene_mappings.tsv gs://data-public-monarchinitiative/monarch-gene-mapping/${RELEASE}/
                    gsutil cp output/gene_mappings.tsv gs://monarch-ingest-data-cache/monarch/
                    
                    gsutil rm gs://data-public-monarchinitiative/monarch-gene-mapping/latest/*
                    gsutil cp gs://data-public-monarchinitiative/monarch-gene-mapping/${RELEASE}/gene_mappings.tsv gs://data-public-monarchinitiative/monarch-gene-mapping/latest/gene_mappings.tsv
                '''
            }
        }
    }
}
