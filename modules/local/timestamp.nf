process TIMESTAMP {
    output:
    env timestamp, emit: timestamp

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    timestamp=\$(date +%s | tr -d '\n\t\r ')
    """
}
