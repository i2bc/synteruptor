FROM nginx:1.23.1-alpine
ADD nginx/default.conf /etc/nginx/conf.d

RUN apk upgrade
RUN apk add git

RUN git clone https://github.com/i2bc/synteruptor_web.git && cd synteruptor_web && git checkout v1.0.1
RUN mkdir -p /var/www/html
RUN mv /synteruptor_web/ /var/www/html
RUN mkdir /var/www/html/synteruptor_web/db
#RUN mv /synteruptor_web /usr/share/nginx/html/

