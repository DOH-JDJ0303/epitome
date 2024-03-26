process SUMMARY {
    label 'process_low'
    container 'docker.io/johnjare/spree:1.0'

    input:
    path clusters
    path lengths
    path ani_ava
    path ani_seeds
    path seeds


    output:
    path "*.csv", emit: summary
    path "*.jpg", emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # remove unwanted headers in cluster dataset
    cat ${clusters} | grep -v 'seq,taxa,segment,cluster' > clusters-no-header.csv
    # run script
    summary.R clusters-no-header.csv ${lengths} ${ani_ava} ${ani_seeds} ${seeds}
    """
}
