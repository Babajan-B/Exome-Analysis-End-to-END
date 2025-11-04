#!/bin/bash
cd $(dirname $0)
java -Xmx4g -jar tools/igv/IGV_2.15.2/igv.jar "$@"
