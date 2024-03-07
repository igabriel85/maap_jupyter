FROM continuumio/miniconda3:4.7.12

RUN mkdir /projects
RUN sed -i -e 's/\/root/\/projects/g' /etc/passwd

RUN conda install -c conda-forge conda=4.7.12 jupyterlab=2.1.4 nodejs=10.13.0 gitpython=3.1.3
RUN conda install -c conda-forge papermill

USER root

# Grant access to jupyterlab config files for base url rewriting
RUN chmod a+rwx -R /opt/conda/lib/python*/site-packages/

# Adjust permissions on /etc/passwd so writable by group root.
RUN chmod g+w /etc/passwd
RUN chmod g+w /etc/environment

###############################
# Custom Jupyter Extensions
###############################
COPY maap-jupyter-ide /maap-jupyter-ide

WORKDIR /maap-jupyter-ide
RUN npm install typescript -g

# control che side panel extension
RUN cd hide_side_panel && npm install \
    && npm run build \
    && jupyter labextension install --no-build

RUN jupyter lab build && \
    jupyter lab clean && \
    jlpm cache clean && \
    npm cache clean --force && \
    rm -rf $HOME/.node-gyp && \
    rm -rf $HOME/.local



#####################################"
# ESA MAAP STACK
##################################

RUN apt-get -y --allow-releaseinfo-change update

RUN  apt-get -y update && \
#      apt-get install -y  --no-install-recommends  rsync software-properties-common python3-lxml python3-pip python3-dev gdal-bin libgdal-dev octave && \
#       apt-get install -y libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 && \
  
      apt-get -y install python3-pip && \
apt-get -y install sudo wget unzip sqlite3 libsqlite3-dev && \
      apt-get -y install jq && \
     apt-get -y install gcc && \
     apt-get -y install g++ && \
     apt-get -y install gdal-bin && \
     apt-get -y install libgdal-dev && \
     pip install requests && \
     pip install OWSLib && \
     apt-get install -y unzip git curl vim && \
     apt -y install gnudatalanguage && \
     apt-get -y install s3fs


#RUN pip install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==`gdal-config --version` && \

# We add the script folder and the zip file to be able to unzip the structure of the project
COPY initTemplate.sh /usr/bmap/initTemplate.sh
COPY initCredentials.sh /usr/bmap/initCredentials.sh
COPY Project_template.zip /usr/bmap/Project_template.zip
COPY shareAlgorithm.sh /usr/bmap/shareAlgorithm.sh
COPY .gitlab-ci.yml /usr/bmap/.gitlab-ci.yml
COPY .condarc /usr/bmap/.condarc

# We add the RestClient file
COPY RestClient.py /usr/bmap/RestClient.py
COPY quicklook_raster.py /usr/bmap/quicklook_raster.py
COPY ingestData.py /usr/bmap/ingestData.py
COPY ingestData.sh /usr/bmap/ingestData.sh
COPY installLib.sh /usr/bmap/.installLib.sh
COPY IN_biosar2.zip /usr/bmap/IN_biosar2.zip
COPY Scripts.zip /usr/bmap/Scripts.zip
RUN chmod a+rwx -R /usr/bmap/

RUN  chmod +x /usr/bmap/initTemplate.sh
RUN  chmod +x /usr/bmap/shareAlgorithm.sh
RUN  chmod +x /usr/bmap/ingestData.sh

ENV PATH="/projects/.maap/bin:/usr/bmap/:${PATH}"
ENV PYTHONPATH="/usr/bmap/:${PYTHONPATH}"
ENV PATH="/usr/modules/:${PATH}"
ENV PYTHONPATH="/usr/modules/:${PYTHONPATH}"

#################################
# END OF ESA CUSTOMISATION
#################################

RUN mkdir /.jupyter
RUN  conda config --set auto_activate_base false && \
    echo "echo 'wait, we create conda environment...'" >> /root/.bashrc && \
    touch /root/.bashrc && \
    echo "cp /usr/bmap/.condarc  /usr/bmap/.installLib.sh /projects" >> /root/.bashrc && \
    echo "conda deactivate" >> /root/.bashrc && \
    echo "if [[ ! -d "/projects/.conda" ]]; then" >>  /root/.bashrc && \
    echo "  conda create -n maap python=3.7.7 -y 2>&1 >/dev/null" >> /root/.bashrc && \
    echo "fi" >> /root/.bashrc && \    
    echo "conda activate maap 2>&1 >/dev/null" >> /root/.bashrc && \

    echo "if [[ ! -d "/projects/.conda" ]]; then" >>  /root/.bashrc && \
    echo "  conda uninstall gdal -y 2>&1 >/dev/null" >> /root/.bashrc && \
    echo "fi" >> /root/.bashrc && \

    echo "conda install ipython_genutils  ipykernel=6.2.0 -y 2>&1 >/dev/null" >> /root/.bashrc && \
    echo "ipython kernel install --user --name=Sar_Training 2>&1 >/dev/null" >> /root/.bashrc && \

    echo "conda install pip -y 2>&1 >/dev/null" >> /root/.bashrc && \
        
    echo "pip install numpy==1.19.1 scipy==1.5.2 requests matplotlib==3.3.1 2>&1 >/dev/null" >> /root/.bashrc && \
    echo "conda install gdal=2.3.3 -y 2>&1 >/dev/null" >> /root/.bashrc && \
    echo "if [[ ! -d "IN_biosar2" ]]; then" >>  /root/.bashrc && \
    echo "unzip /usr/bmap/IN_biosar2.zip  2>&1 >/dev/null" >> /root/.bashrc && \
    echo "fi" >> /root/.bashrc && \

    echo "if [[ ! -d "Scripts" ]]; then" >>  /root/.bashrc && \
    echo "unzip /usr/bmap/Scripts.zip -d /projects/Scripts  2>&1 >/dev/null" >> /root/.bashrc && \
    echo "fi" >> /root/.bashrc && \

    echo "./.installLib.sh" >> /root/.bashrc && \
    echo "echo 'Workspace is ready'" >> /root/.bashrc && \
    echo "cd /projects >& /dev/null" >> /root/.bashrc && \
    echo "export PATH=".conda/envs/maap/bin:${PATH}"" >> /root/.bashrc && \    

    echo "if [[ ! -f ".init" ]]; then" >>  /root/.bashrc && \
    echo "touch .init" >> /root/.bashrc && \
    echo "exit" >> /root/.bashrc && \
    echo "fi" >> /root/.bashrc

RUN export PATH="/usr/bmap/:$PATH"
RUN echo "PATH=/opt/conda/bin:${PATH}:/usr/bmap/" >> /etc/environment

RUN chmod a+rwx -R /maap-jupyter-ide/ && chmod a+rwx -R /.jupyter

ENV SHELL="bash"

WORKDIR /projects
EXPOSE 3100

ENTRYPOINT ["/bin/bash", "/maap-jupyter-ide/entrypoint.sh"]
