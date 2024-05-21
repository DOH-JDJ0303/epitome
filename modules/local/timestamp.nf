process TIMESTAMP {
    container = 'docker.io/jdj0303/epitome-base:1.0.0'

    output:
    env timestamp, emit: timestamp

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    timestamp=\$(date +%s | tr -d '\n\t\r ')
    """
}
