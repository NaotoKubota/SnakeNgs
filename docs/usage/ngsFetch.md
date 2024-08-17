# ngsFetch

Fetching NGS data from the public server using [ffq](https://github.com/pachterlab/ffq) and [aria2](https://aria2.github.io/).

## Usage

``` bash
$ docker run --rm naotokubota/ngsFetch:latest ngsFetch -h

                       ____        __           __         
                      /\  _`\     /\ \__       /\ \        
  ___      __     ____\ \ \L\_\ __\ \ ,_\   ___\ \ \___    
/' _ `\  /'_ `\  /',__\\ \  _\/'__`\ \ \/  /'___\ \  _ `\  
/\ \/\ \/\ \L\ \/\__, `\\ \ \/\  __/\ \ \_/\ \__/\ \ \ \ \ 
\ \_\ \_\ \____ \/\____/ \ \_\ \____\\ \__\ \____\\ \_\ \_\
 \/_/\/_/\/___L\ \/___/   \/_/\/____/ \/__/\/____/ \/_/\/_/
           /\____/                                         
           \_/__/                                          

Version: v0.1


ngsFetch v0.1 -Pipeline for fetching NGS data in ZhengLab-

Usage: ngsFetch -i ID -o output_dir

    -h  Display help
    -v  Show version
    -i  ID
    -o  Output directory
```

You are supposed to get the data from the public server using the ID provided by the user. The data will be stored in the output directory after checking the [MD5](https://en.wikipedia.org/wiki/MD5) to ensure the integrity of the data.

Please refer to the [tutorial](../tutorial/Fetching.md) for more information.

## Software

- [ffq](https://github.com/pachterlab/ffq) (v0.3.1)
- [aria2](https://aria2.github.io/)

## Docker image

- [naotokubota/ngsFetch](https://hub.docker.com/r/naotokubota/ngsFetch)
