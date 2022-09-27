pipeline {
    agent { label 'monarch-agent-small' }
        environment {
        HOME = "${env.WORKSPACE}"
        RELEASE = sh(script: "echo `date +%Y-%m-%d`", returnStdout: true).trim()
        PATH = "/home/ubuntu/.poetry/bin:${env.PATH}"
    }
    stages {
        stage('setup') {
            steps {
                sh 'poetry install'
            }
        }
        stage('generate-mapping-file') {
            steps {
                sh 'poetry run gene-mapping generate --download'
            }
        }
        stage('upload-mapping-file'){
            steps{
                sh 'gsutil cp output/gene_mappings.tsv gs://monarch-archive/monarch-gene-mapping/${RELEASE}/'
            }
        }
    }
}
