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
                    gsutil cp output/gene_mappings.sssom.tsv gs://monarch-archive/monarch-gene-mapping/${RELEASE}/
                    gsutil cp output/gene_mappings.sssom.tsv gs://data-public-monarchinitiative/monarch-gene-mapping/${RELEASE}/
                    gsutil cp output/gene_mappings.sssom.tsv gs://monarch-ingest-data-cache/monarch/
                    
                    gsutil rm gs://data-public-monarchinitiative/monarch-gene-mapping/latest/*
                    gsutil cp gs://data-public-monarchinitiative/monarch-gene-mapping/${RELEASE}/gene_mappings.sssom.tsv gs://data-public-monarchinitiative/monarch-gene-mapping/latest/gene_mappings.sssom.tsv
                '''
            }
        }
        stage('index') {
            steps {
                sh '''
                    echo "Current directory: $(pwd)"
                    python3 --version
                    pip --version
                    export PATH=$HOME/.local/bin:$PATH
                    echo "Path: $PATH"

                    cd $HOME
                    mkdir data-public
                    gcsfuse --implicit-dirs data-public-monarchinitiative data-public

                    git clone https://github.com/monarch-initiative/monarch-file-server.git
                    pip install -r monarch-file-server/scripts/requirements.txt
                    python3 monarch-file-server/scripts/directory_indexer.py --inject monarch-file-server/scripts/directory-index-template.html --directory data-public --prefix https://data.monarchinitiative.org -x
                '''
            }
        }
    }
}
