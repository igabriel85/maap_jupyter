`docker build -t <registryurl>/root/jupyter_image/gchang-maap-che7-esa  .`

`docker push <registryurl>/root/jupyter_image/gchang-maap-che7-esa:latest`

```
apiVersion: 1.0.0
metadata:
  generateName: develop-
attributes:
  editorFree: 'true'
components:
  - endpoints:
      - attributes:
          type: ide
          discoverable: 'false'
          path: /
          protocol: http
          public: 'true'
        name: jupyter
        port: 3100
    referenceContent: |
      kind: List
      items:
        - apiVersion: v1
          kind: Pod
          metadata:
            name: ws
            labels:
              name: ws
          spec:
           containers: 
            - name: jupyter
              image: '<registryurl>/root/jupyter_image/gchang-maap-che7-esa:latest'
              imagePullPolicy: Always
              resources:
                limits:
                  memory: 2048Mi
              securityContext:
                privileged: true
    type: kubernetes
    alias: maap-jupyterlab
    env:
      - value: WS_JUPYTER
        name: MACHINE_NAME

```
