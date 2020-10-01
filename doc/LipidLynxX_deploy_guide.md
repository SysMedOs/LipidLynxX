# Deploy LipidLynxX

## Requirement

LipidLynxX is based on FastAPI, thus tool chains using `nginx`, `supervisor`,`gunicorn`, and `uvicorn` is recommended.

+ **System installation**
  
  + `sudo apt install nginx supervisor`

+ **Python installation**
  
  + Install `anaconda`/`miniconda` and create a virtual enviroment e.g. `envlynx`
    
    + `conda create -n envlynx python=3.7`
    
    + `conda activate envlynx`
    
    + `pip install -r requirement.txt`  # do this under LipidLynxX folder
  
  + Install `gunicorn`, and `uvicorn`
    
    + `pip install gunicorn uvicorn gevent`



Restart linux server before going to following setps.



## Setup nginx

- Check if you have `nginx` installed and go to its config folder
  
  - `cd /etc/nginx/conf.d`

- Create a new config file named `lynx.conf`
  
  - `sudo nano lynx.conf`

- Edit following content accordingly especially content inside`{}` and save the file.
  
  - ```json
    server {
        listen 80;
        server_name {example.com};
        access_log  /var/log/nginx/example.log;
    
        location / {
            proxy_pass http://127.0.0.1:8000;
            proxy_set_header Host $host;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        }
    }
    
    ```
    
    - `http://127.0.0.1:8000` is default by `gunicorn`
- Restart `nginx`
  - `sudo service nginx restart`

## Setup supervisor

+ Check if you have `supervisor` installed and go to its config folder
  
  + `cd /etc/supervisor/conf.d`

+ Create a new config file named `lynx.conf`
  
  + `sudo nano lynx.conf`

+ Edit following content accordingly especially content inside`{}` and save the file.
  
  + ```ini
    [program:daemon_lynx]
    directory=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}
    environment=PATH=/home/anaconda3/envs/envlynx38/bin  # or your env path
    command=/home/{USER_NAME}/anaconda3/envs/envlynx38/bin/python daemon_lynx.py
    autorestart=true
    redirect_stderr=true
    stdout_logfile=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}/daemon_lynx_stdout.log
    stderr_logfile=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}/daemon_lynx_error.log
    
    [program:lynx_app]
    directory=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}
    environment=PATH=/home/anaconda3/envs/envlynx38/bin  # or your env path
    command=/home/{USER_NAME}/anaconda3/envs/envlynx38/bin/gunicorn lynx.app:app -w 4 -k uvicorn.workers.UvicornWorker
    autorestart=true
    redirect_stderr=true
    stdout_logfile=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}/lynx_app_stdout.log
    stderr_logfile=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}/lynx_app_error.log

    [supervisord]
    logfile=/home/{USER_NAME}/{PATH_TO_LipidLynxX_FOLDER}/supervisord.log
    logfile_maxbytes=50MB

    ```
    
    + if you use miniconda, then replace anaconda3 by miniconda/miniconda3
    
    + `-w 4` means set up 4 workers, you can adjust number of workers based on your server configuration. 

+ Refresh configs for `supervisor`
  
  + `sudo supervisorctl reread`

+ Restart `supervisor`
  
  + `sudo service supervisor restart`



## Test LipidLynxX

Now you can visit `127.0.0.1` on the PC/server running LipidLynxX and other PCs in the same local network can access LipidLynxX web GUI using the IP of LipidLynxX computer e.g. `192.168.100.101`



## Update LipidLynxX

After updating LipidLynxX, you have to restart the web service.

You can reboot the server or run following two commands to restart `supervisor` and `nginx` service.

+ stop supervisor `sudo supervisorctl stop lynx_app lynx_daemon`
+ check supervisor status `sudo supervisorctl status`
+ kill rest of daemon `sudo pkill -f daemon_lynx`
+ start supervisor `sudo supervisorctl start lynx_app lynx_daemon`
+ check supervisor status `sudo supervisorctl status`
+ restart nginx `sudo service nginx restart`

LipidLynxX should then run with the new version.
